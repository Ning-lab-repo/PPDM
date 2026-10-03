# PPDM cross-validation sensitivity analysis

Bolt-on script for the hyperparameter-sensitivity analysis reported in the manuscript.
It is deliberately **standalone**: the training pipeline, the checkpoints and the reported
results are never modified, and **no third-party tuning framework** (Optuna,
`GridSearchCV`, `RandomizedSearchCV`) is used — the cross-validation logic is written out
in full.

## Why it is a separate script

The final models were **not** tuned per endpoint. Hyperparameters were pre-specified
(hidden size 32, batch size 128, learning rate 1 × 10⁻³, 20 epochs) and shared by all 370
disease-specific models. This script characterises how sensitive that choice is; it does
not select hyperparameters.

## Partition discipline

```
England participants
├── 80 % "training part"        ← 3-fold stratified CV runs ONLY inside here
└── 20 % internal validation    ← never touched
Wales + Scotland (external test) ← never touched
```

The 80/20 split reproduces the training code exactly:

```python
train_test_split(idx, test_size=0.2, random_state=42, stratify=y)
```

so the CV estimate is directly comparable with the reported validation performance and
cannot leak into it.

## Files

| File | Purpose |
|---|---|
| `ppdm_cv_sensitivity.py` | The complete analysis: model definition + CV harness + CLI |
| `README.md` | This file |

## Model definition

`LSTM_SelfAttention_Model` in `ppdm_cv_sensitivity.py` is the class from the PPDM
training code, instantiated with `use_time_aware=False, use_dnn=False` (the behaviour-only
variant reported in the paper). Verified against the training code:

- identical parameter names (66) and shapes, 69,889 parameters
- identical forward output on the same input (**maximum absolute difference 0.0**)
- identical attention weights

## Required input files

These are derived from UK Biobank data and therefore **cannot be redistributed**; point
`--data-dir` at a directory containing them.

| File | Shape / columns | Notes |
|---|---|---|
| `movement_imputed.npy` | `(N, 48, 4)` float32 | hourly behavioural-state composition; rows 0–23 = weekday hours, rows 24–47 = weekend hours; columns = SB, LPA, MVPA, sleep; each row sums to 1 |
| `labels_370.npy` | `(N, 370)` float32 | per-endpoint incident-disease labels (0 / 1); NaN where not evaluable |
| `labels_meta.csv` | `disease, j, n_case_england, n_case_test` | maps an ICD-10 level-3 code to its column `j` |
| `movement_ids.csv` | `Participant ID, test` | `test == 1` → England (training partition); `test == 2` or `3` → Wales / Scotland (external test) |

`N` must be identical across all four files.

## Usage

```bash
# smoke test on synthetic data — no input files needed
python ppdm_cv_sensitivity.py --self-test

# one run: one endpoint, one configuration, one fold
python ppdm_cv_sensitivity.py --data-dir ./data --outdir ./hp_runs \
    --disease I10 --hidden 32 --batch 128 --lr 1e-3 --fold 0

# all 27 configurations × 3 folds for one endpoint
python ppdm_cv_sensitivity.py --data-dir ./data --outdir ./hp_runs --disease I10 --grid

# the full sensitivity analysis: 20 endpoints × 27 configurations × 3 folds = 1,620 runs
python ppdm_cv_sensitivity.py --data-dir ./data --outdir ./hp_runs --all

# list the 20 representative endpoints
python ppdm_cv_sensitivity.py --list-endpoints
```

Runs are **resumable**: a run whose output file already exists is skipped, so an
interrupted sweep can simply be restarted. Each run is single-threaded
(`torch.set_num_threads(1)`), so the number of concurrent processes equals the number of
cores used.

### Scale and runtime

The full sweep is 20 endpoints × 27 configurations × 3 folds = **1,620 runs**. A single
run trains the model once on roughly 80 % of the England training part (≈ 37,000
participants for a common endpoint) for 20 epochs, which takes on the order of 10–25
minutes per run on CPU depending on endpoint prevalence. Because every run is
single-threaded the sweep is embarrassingly parallel: the reported analysis was executed
with 60–120 concurrent processes and completed in about 5 hours. For a quick check, use
`--self-test`, which needs no input files.

## Output

One JSON per run: `<outdir>/<disease>_h<hidden>_b<batch>_lr<lr>_f<fold>.json`

```json
{"disease": "I10", "hidden": 32, "batch": 128, "lr": 0.001, "fold": 0, "folds": 3,
 "epochs": 20, "seed": 42, "n_fit": 37168, "n_val": 18584, "n_case_val": 2661,
 "AUC": 0.6079, "AUPRC": 0.2002, "prev_val_pct": 14.3188, "seconds": 1550.5}
```

## Grid and representative endpoints

| Dimension | Values | Pre-specified value |
|---|---|---|
| hidden size | 16, **32**, 64 | 32 (centre of the grid) |
| batch size | 64, **128**, 256 | 128 (centre) |
| learning rate | 5 × 10⁻⁴, **1 × 10⁻³**, 2 × 10⁻³ | 1 × 10⁻³ (centre) |

All three pre-specified values are the **centre** of the grid, so every conclusion about
"the pre-specified configuration lies within the observed range" is supported on both
sides.

The 20 representative endpoints (14 ICD chapters, prevalence 0.126 %–14.32 %) were
**selected before any hyperparameter evaluation**: ten endpoints highlighted in the paper,
plus ten drawn by a fixed rule (prevalence tertile × ICD-chapter diversity).

## Environment

`torch`, `numpy`, `pandas`, `scikit-learn`. CPU only; no GPU required.
