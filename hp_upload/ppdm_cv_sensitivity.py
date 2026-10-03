"""PPDM cross-validation sensitivity analysis -- self-contained, bolt-on script.

This script answers one question only: for a given endpoint and hyperparameter
configuration, how well does the PPDM architecture discriminate under cross-validation
inside the England training partition?  It is deliberately standalone: the training
pipeline, the checkpoints and the reported results are never touched, and no
third-party tuning framework (Optuna / GridSearchCV / RandomizedSearchCV) is used.

----------------------------------------------------------------------------------------
Partition discipline
----------------------------------------------------------------------------------------
    England participants
        |-- 80 % "training part"        <-- k-fold stratified CV runs ONLY inside here
        |-- 20 % internal validation    <-- never touched
    Wales + Scotland (external test)    <-- never touched

The 80/20 split reproduces the training code exactly
(`train_test_split(..., test_size=0.2, random_state=42, stratify=y)`), so the CV
estimate is directly comparable with the reported validation performance and cannot
leak into it.

----------------------------------------------------------------------------------------
Required input files
----------------------------------------------------------------------------------------
These are derived from UK Biobank data and therefore cannot be redistributed; point
`--data-dir` at a directory containing them.

    movement_imputed.npy   float32, shape (N, 48, 4)
                           hourly behavioural-state composition: the first 24 rows are the
                           weekday hours and the last 24 rows are the weekend hours; the 4
                           columns are SB, LPA, MVPA and sleep.  Rows sum to 1 within each hour.
    labels_370.npy         float32, shape (N, 370)
                           per-endpoint incident-disease labels (0 / 1); NaN where the
                           participant is not evaluable for that endpoint.
    labels_meta.csv        columns: disease, j, n_case_england, n_case_test
                           maps an ICD-10 level-3 code to its column index j in labels_370.npy
    movement_ids.csv       columns: Participant ID, test
                           test == 1 -> England (training partition); test == 2 or 3 ->
                           Wales / Scotland (external test partition)

N must be identical across the two .npy files and the two .csv files.

----------------------------------------------------------------------------------------
Usage
----------------------------------------------------------------------------------------
    # one run (one endpoint, one configuration, one fold)
    python ppdm_cv_sensitivity.py --data-dir ./data --outdir ./hp_runs \
        --disease I10 --hidden 32 --batch 128 --lr 1e-3 --fold 0

    # all 27 configurations x 3 folds for one endpoint
    python ppdm_cv_sensitivity.py --data-dir ./data --outdir ./hp_runs \
        --disease I10 --grid

    # the full sensitivity analysis: 20 representative endpoints x 27 configurations x 3 folds
    python ppdm_cv_sensitivity.py --data-dir ./data --outdir ./hp_runs --all

    # smoke test on synthetic data (no UK Biobank access needed)
    python ppdm_cv_sensitivity.py --self-test

Each run writes one JSON file and can be re-run safely: runs whose output file already
exists are skipped, so an interrupted sweep simply resumes.

----------------------------------------------------------------------------------------
Output
----------------------------------------------------------------------------------------
    <outdir>/<disease>_h<hidden>_b<batch>_lr<lr>_f<fold>.json
    fields: disease, hidden, batch, lr, fold, folds, epochs, seed,
            n_fit, n_val, n_case_val, AUC, AUPRC, prev_val_pct, seconds

Environment: torch, numpy, pandas, scikit-learn.  CPU only; each run is single-threaded.
"""
import argparse
import itertools
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold, train_test_split
from torch.utils.data import DataLoader, TensorDataset

# --------------------------------------------------------------------------------------
# Model definition.
# Copied verbatim from the training code (`LSTM_SelfAttention_Model` in the PPDM training
# notebooks); the cross-validation harness below instantiates it with
# use_time_aware=False and use_dnn=False, i.e. the behaviour-only variant reported in the
# paper, so the architecture is identical to the one used for the published models.
# --------------------------------------------------------------------------------------
class LSTM_SelfAttention_Model(nn.Module):
    def __init__(self, input_size, hidden_size, output_size, static_feature_size,
                 lstm_dropout_rate=0.5, attention_dropout_rate=0.3,
                 use_time_aware=True, use_dnn=True):
        super(LSTM_SelfAttention_Model, self).__init__()
        self.use_time_aware = use_time_aware
        self.use_dnn = use_dnn

        self.lstm_workday = nn.ModuleList([nn.LSTM(input_size=1, hidden_size=hidden_size, batch_first=True)
                                           for _ in range(input_size)])
        self.lstm_restday = nn.ModuleList([nn.LSTM(input_size=1, hidden_size=hidden_size, batch_first=True)
                                           for _ in range(input_size)])

        self.lstm_dropout = nn.Dropout(lstm_dropout_rate)
        self.attention_dropout = nn.Dropout(attention_dropout_rate)

        self.self_attention_workday = nn.ModuleList([nn.MultiheadAttention(embed_dim=hidden_size, num_heads=2)
                                                     for _ in range(input_size)])
        self.self_attention_restday = nn.ModuleList([nn.MultiheadAttention(embed_dim=hidden_size, num_heads=2)
                                                     for _ in range(input_size)])

        if self.use_time_aware:
            self.time_aware_module = nn.Sequential(
                nn.Linear(1, 16), nn.ReLU(), nn.Dropout(0.5),
                nn.Linear(16, 8), nn.ReLU(), nn.Dropout(0.5),
            )

        if self.use_dnn:
            self.dnn = nn.Sequential(
                nn.Linear(static_feature_size, 64), nn.ReLU(), nn.Dropout(0.5),
                nn.Linear(64, 32), nn.ReLU(), nn.Dropout(0.5),
            )

        fc_input_size = hidden_size * 2 * input_size
        if self.use_dnn:
            fc_input_size += 32
        if self.use_time_aware:
            fc_input_size += 8
        self.fc = nn.Linear(fc_input_size, output_size)

    def self_attention(self, lstm_out, self_attention_layer):
        lstm_out = lstm_out.permute(1, 0, 2)
        attn_output, attn_weights = self_attention_layer(lstm_out, lstm_out, lstm_out)
        attn_output = attn_output.permute(1, 0, 2)
        return torch.mean(attn_output, dim=1), attn_weights

    def forward(self, workday_input, restday_input, static_input, time_gap):
        workday_attended_list, workday_attn_weights_list = [], []
        for i in range(workday_input.size(2)):
            out, _ = self.lstm_workday[i](workday_input[:, :, i].unsqueeze(2))
            out = self.lstm_dropout(out)
            attended, w = self.self_attention(out, self.self_attention_workday[i])
            workday_attended_list.append(self.attention_dropout(attended))
            workday_attn_weights_list.append(w)

        restday_attended_list, restday_attn_weights_list = [], []
        for i in range(restday_input.size(2)):
            out, _ = self.lstm_restday[i](restday_input[:, :, i].unsqueeze(2))
            out = self.lstm_dropout(out)
            attended, w = self.self_attention(out, self.self_attention_restday[i])
            restday_attended_list.append(self.attention_dropout(attended))
            restday_attn_weights_list.append(w)

        combined = torch.cat(workday_attended_list + restday_attended_list, dim=1)
        combined = self.lstm_dropout(combined)

        if self.use_dnn:
            combined = torch.cat((combined, self.dnn(static_input)), dim=1)
        if self.use_time_aware:
            combined = torch.cat((combined, self.time_aware_module(time_gap)), dim=1)

        return self.fc(combined), workday_attn_weights_list, restday_attn_weights_list


# --------------------------------------------------------------------------------------
# Representative endpoints used for the sensitivity analysis.
# Selected before any hyperparameter evaluation: ten endpoints highlighted in the paper,
# plus ten drawn by a fixed rule (prevalence tertile x ICD-chapter diversity).
# --------------------------------------------------------------------------------------
DEFAULT_ENDPOINTS = [
    "I10", "K57", "E66", "F32", "E11", "N18", "N17", "J44", "B95", "D70",
    "L40", "H43", "G62", "G30", "C64", "J38", "C85", "N47", "M40", "H83",
]

# 3 x 3 x 3 candidate configurations; the pre-specified configuration
# (hidden 32, batch 128, lr 1e-3) is the centre of the grid.
GRID_HIDDEN = [16, 32, 64]
GRID_BATCH = [64, 128, 256]
GRID_LR = [5e-4, 1e-3, 2e-3]
PRESPECIFIED = (32, 128, 1e-3)


# --------------------------------------------------------------------------------------
# Data
# --------------------------------------------------------------------------------------
class Data:
    """Loads the cached behaviour matrices and labels and rebuilds the training split."""

    def __init__(self, data_dir):
        self.dir = data_dir

        def need(name):
            p = os.path.join(data_dir, name)
            if not os.path.exists(p):
                raise SystemExit(f"missing input file: {p}")
            return p

        self.X = np.load(need("movement_imputed.npy"))          # (N, 48, 4)
        self.labels = np.load(need("labels_370.npy"))           # (N, 370)
        self.meta = pd.read_csv(need("labels_meta.csv"))
        ids = pd.read_csv(need("movement_ids.csv"))

        if self.X.shape[0] != self.labels.shape[0] or self.X.shape[0] != len(ids):
            raise SystemExit("input files disagree on the number of participants")
        if self.X.shape[1] != 48 or self.X.shape[2] != 4:
            raise SystemExit(f"movement_imputed.npy must be (N, 48, 4), got {self.X.shape}")

        # test == 1 -> England (training partition); test == 2/3 -> Wales/Scotland (external)
        self.england = np.where(ids["test"].to_numpy() == 1)[0]
        self.col = dict(zip(self.meta.disease.astype(str), self.meta.j.astype(int)))

    def endpoint(self, disease):
        if disease not in self.col:
            raise SystemExit(f"unknown endpoint {disease!r}; see labels_meta.csv")
        y_all = self.labels[:, self.col[disease]]
        keep = np.isin(y_all[self.england], [0, 1])
        eng_ok = self.england[keep]
        return eng_ok, y_all[eng_ok].astype(np.float32)


# --------------------------------------------------------------------------------------
# One CV run
# --------------------------------------------------------------------------------------
def run_one(data, disease, hidden, batch, lr, fold, folds, epochs, seed, device="cpu"):
    torch.set_num_threads(1)
    torch.manual_seed(seed)
    np.random.seed(seed)
    t0 = time.time()

    eng_ok, y = data.endpoint(disease)

    # reproduce the training code's 8:2 split and keep only the 80 % training part
    idx = np.arange(len(y))
    tr_idx, _va_idx = train_test_split(idx, test_size=0.2, random_state=42, stratify=y)

    skf = StratifiedKFold(n_splits=folds, shuffle=True, random_state=42)
    fit_pos, val_pos = list(skf.split(tr_idx, y[tr_idx]))[fold]
    fit_idx, val_idx = tr_idx[fit_pos], tr_idx[val_pos]

    Xw = data.X[eng_ok][:, :24, :]
    Xr = data.X[eng_ok][:, 24:, :]

    wd_t = torch.tensor(np.ascontiguousarray(Xw[fit_idx]), dtype=torch.float32)
    rd_t = torch.tensor(np.ascontiguousarray(Xr[fit_idx]), dtype=torch.float32)
    y_t = torch.tensor(y[fit_idx]).unsqueeze(1)
    loader = DataLoader(TensorDataset(wd_t, rd_t, y_t), batch_size=batch, shuffle=True)

    # behaviour-only variant: no static-feature encoder, no time-gap encoder
    model = LSTM_SelfAttention_Model(input_size=4, hidden_size=hidden, output_size=1,
                                     static_feature_size=1, use_time_aware=False, use_dnn=False)
    model.to(device)
    criterion = nn.BCEWithLogitsLoss()
    optimizer = optim.Adam(model.parameters(), lr=lr)

    model.train()
    for _ in range(epochs):
        for w, r, yy in loader:
            optimizer.zero_grad()
            outputs, _, _ = model(w, r, None, None)
            criterion(outputs, yy).backward()
            optimizer.step()

    model.eval()
    preds = []
    with torch.no_grad():
        for i in range(0, len(val_idx), 8192):
            w = torch.tensor(np.ascontiguousarray(Xw[val_idx][i:i + 8192]), dtype=torch.float32)
            r = torch.tensor(np.ascontiguousarray(Xr[val_idx][i:i + 8192]), dtype=torch.float32)
            outputs, _, _ = model(w, r, None, None)
            preds.append(torch.sigmoid(outputs).squeeze(1).cpu().numpy())
    p = np.concatenate(preds)
    yv = y[val_idx].astype(int)

    return {
        "disease": disease, "hidden": hidden, "batch": batch, "lr": lr,
        "fold": fold, "folds": folds, "epochs": epochs, "seed": seed,
        "n_fit": int(len(fit_idx)), "n_val": int(len(val_idx)), "n_case_val": int(yv.sum()),
        "AUC": float(roc_auc_score(yv, p)), "AUPRC": float(average_precision_score(yv, p)),
        "prev_val_pct": round(100 * float(yv.mean()), 4),
        "seconds": round(time.time() - t0, 1),
    }


def tag_of(disease, hidden, batch, lr, fold):
    return f"{disease}_h{hidden}_b{batch}_lr{lr:g}_f{fold}"


def run_and_save(data, outdir, **kw):
    os.makedirs(outdir, exist_ok=True)
    tag = tag_of(kw["disease"], kw["hidden"], kw["batch"], kw["lr"], kw["fold"])
    path = os.path.join(outdir, tag + ".json")
    if os.path.exists(path):
        print(f"[skip] {tag}", flush=True)
        return
    rec = run_one(data, **kw)
    with open(path, "w") as fh:
        json.dump(rec, fh)
    print(json.dumps(rec), flush=True)


# --------------------------------------------------------------------------------------
# Self test on synthetic data
# --------------------------------------------------------------------------------------
def self_test(n=1200, seed=0, outdir=None):
    """Generates a small synthetic cache and runs one fold end to end."""
    outdir = outdir or os.path.join(os.getcwd(), "_self_test")
    ddir = os.path.join(outdir, "data")
    rdir = os.path.join(outdir, "hp_runs")
    os.makedirs(ddir, exist_ok=True)
    rng = np.random.default_rng(seed)

    X = rng.random((n, 48, 4)).astype(np.float32)
    X /= X.sum(axis=2, keepdims=True)
    y = (rng.random((n, 370)) < 0.05).astype(np.float32)
    np.save(os.path.join(ddir, "movement_imputed.npy"), X)
    np.save(os.path.join(ddir, "labels_370.npy"), y)
    pd.DataFrame({"disease": ["SYN01", "SYN02"], "j": [0, 1],
                  "n_case_england": [int(y[:, 0].sum()), int(y[:, 1].sum())],
                  "n_case_test": [0, 0]}).to_csv(os.path.join(ddir, "labels_meta.csv"), index=False)
    pd.DataFrame({"Participant ID": np.arange(n),
                  "test": np.where(np.arange(n) % 10 < 9, 1, 2)}).to_csv(
        os.path.join(ddir, "movement_ids.csv"), index=False)

    data = Data(ddir)
    rec = run_one(data, "SYN01", hidden=8, batch=128, lr=1e-3, fold=0,
                  folds=3, epochs=2, seed=42)
    print("self-test OK")
    print(json.dumps(rec, indent=2))
    print(f"\nsynthetic data written to {ddir}")
    return rec


# --------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="PPDM cross-validation sensitivity analysis (bolt-on; the training pipeline is not modified).",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", default="./data",
                    help="directory holding movement_imputed.npy, labels_370.npy, "
                         "labels_meta.csv and movement_ids.csv")
    ap.add_argument("--outdir", default="./hp_runs", help="directory for the per-run JSON files")
    ap.add_argument("--disease", help="ICD-10 level-3 endpoint code, e.g. I10")
    ap.add_argument("--hidden", type=int, help="LSTM hidden size")
    ap.add_argument("--batch", type=int, help="mini-batch size")
    ap.add_argument("--lr", type=float, help="Adam learning rate")
    ap.add_argument("--fold", type=int, help="fold index to evaluate")
    ap.add_argument("--folds", type=int, default=3, help="number of CV folds (default 3)")
    ap.add_argument("--epochs", type=int, default=20, help="training epochs per fold (default 20)")
    ap.add_argument("--seed", type=int, default=42, help="random seed (default 42)")
    ap.add_argument("--grid", action="store_true",
                    help="run all 27 configurations x --folds for --disease")
    ap.add_argument("--all", action="store_true",
                    help="run the full sweep over the 20 representative endpoints")
    ap.add_argument("--endpoints", nargs="*", default=None,
                    help="override the endpoint list used by --all")
    ap.add_argument("--list-endpoints", action="store_true", help="print the endpoint list and exit")
    ap.add_argument("--self-test", action="store_true",
                    help="run a smoke test on synthetic data (no input files required)")
    a = ap.parse_args()

    if a.list_endpoints:
        print(f"{len(DEFAULT_ENDPOINTS)} representative endpoints:")
        print(" ".join(DEFAULT_ENDPOINTS))
        return

    if a.self_test:
        self_test(outdir=a.outdir if a.outdir != "./hp_runs" else None)
        return

    data = Data(a.data_dir)

    if a.all:
        endpoints = a.endpoints or DEFAULT_ENDPOINTS
        for disease in endpoints:
            for hidden, batch, lr in itertools.product(GRID_HIDDEN, GRID_BATCH, GRID_LR):
                for fold in range(a.folds):
                    run_and_save(data, a.outdir, disease=disease, hidden=hidden, batch=batch,
                                 lr=lr, fold=fold, folds=a.folds, epochs=a.epochs, seed=a.seed)
        return

    if a.disease is None:
        ap.error("--disease is required unless --all / --self-test / --list-endpoints is used")

    if a.grid:
        for hidden, batch, lr in itertools.product(GRID_HIDDEN, GRID_BATCH, GRID_LR):
            for fold in range(a.folds):
                run_and_save(data, a.outdir, disease=a.disease, hidden=hidden, batch=batch,
                             lr=lr, fold=fold, folds=a.folds, epochs=a.epochs, seed=a.seed)
        return

    missing = [n for n in ("hidden", "batch", "lr", "fold") if getattr(a, n) is None]
    if missing:
        ap.error("missing required option(s): " + ", ".join("--" + m for m in missing))
    run_and_save(data, a.outdir, disease=a.disease, hidden=a.hidden, batch=a.batch,
                 lr=a.lr, fold=a.fold, folds=a.folds, epochs=a.epochs, seed=a.seed)


if __name__ == "__main__":
    sys.exit(main())
