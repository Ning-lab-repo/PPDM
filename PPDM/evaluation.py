"""
Model evaluation and prediction module
Handles batch prediction, threshold optimization, and performance metrics
"""

import torch
import numpy as np
import pandas as pd
import joblib
import os
from sklearn.metrics import (
    roc_auc_score, precision_recall_curve, auc, 
    f1_score, roc_curve, confusion_matrix
)

class ModelEvaluator:
    def __init__(self, device='cuda' if torch.cuda.is_available() else 'cpu'):
        """Initialize ModelEvaluator with device settings"""
        self.device = torch.device(device)
        
    def load_models(self, model_dir):
        """Load all trained models from directory"""
        models = {}
        if not os.path.exists(model_dir):
            print(f"Model directory not found: {model_dir}")
            return models
            
        for filename in os.listdir(model_dir):
            if filename.endswith(".joblib") or filename.endswith(".pth"):
                disease = filename.replace("_model.joblib", "").replace("_model.pth", "")
                model_path = os.path.join(model_dir, filename)
                
                try:
                    models[disease] = joblib.load(model_path)
                    print(f"Loaded model for {disease}")
                except:
                    print(f"Failed to load model for {disease}")
                    
        return models
    
    def predict_in_batches(self, models, movement_data, participant_ids=None, batch_size=5000):
        """Generate predictions in batches to handle large datasets"""
        n_samples = len(movement_data)
        all_results = []
        
        for i in range(0, n_samples, batch_size):
            batch_end = min(i + batch_size, n_samples)
            batch_data = movement_data[i:batch_end]
            
            # Prepare tensors
            X_tensor = torch.tensor(batch_data, dtype=torch.float32).to(self.device)
            workday_tensor = X_tensor[:, :24, :]
            restday_tensor = X_tensor[:, 24:, :]
            
            # Create dummy tensors for static features and time gap
            # In real usage, replace with actual data
            batch_size_current = len(batch_data)
            static_tensor = torch.zeros((batch_size_current, 1), dtype=torch.float32).to(self.device)
            time_gap_tensor = torch.zeros((batch_size_current, 1), dtype=torch.float32).to(self.device)
            
            # Generate predictions for each disease
            batch_results = {}
            for disease, model in models.items():
                model = model.to(self.device)
                model.eval()
                with torch.no_grad():
                    y_pred_proba, _, _ = model(workday_tensor, restday_tensor, 
                                               static_tensor, time_gap_tensor)
                    y_pred_proba = torch.sigmoid(y_pred_proba.squeeze())
                batch_results[disease] = y_pred_proba.cpu().numpy()
            
            # Create DataFrame for batch results
            batch_df = pd.DataFrame(batch_results)
            if participant_ids is not None:
                batch_df.insert(0, 'Participant ID', participant_ids.iloc[i:batch_end].values)
            
            all_results.append(batch_df)
            print(f"Processed batch {i//batch_size + 1}/{(n_samples-1)//batch_size + 1}")
        
        # Combine all batch results
        final_df = pd.concat(all_results, ignore_index=True)
        return final_df
    
    def find_optimal_threshold(self, y_true, y_pred_proba):
        """Find optimal threshold that maximizes F1 score"""
        precision, recall, thresholds = precision_recall_curve(y_true, y_pred_proba)
        f1_scores = 2 * (precision * recall) / (precision + recall + 1e-8)
        f1_scores = np.nan_to_num(f1_scores)
        
        optimal_idx = np.argmax(f1_scores)
        optimal_threshold = thresholds[optimal_idx] if optimal_idx < len(thresholds) else 0.5
        return optimal_threshold, f1_scores[optimal_idx]
    
    def evaluate_all_metrics(self, risk_predictions, true_labels, output_dir=None):
        """Comprehensive evaluation with multiple metrics"""
        # Reset indices for alignment
        risk_predictions = risk_predictions.reset_index(drop=True)
        true_labels = true_labels.reset_index(drop=True)
        
        # Create output directory if specified
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
        
        # Storage for results
        auc_results = []
        roc_curve_data = []
        pr_curve_data = []
        threshold_results = []
        
        for disease in risk_predictions.columns:
            if disease == "Participant ID":
                continue
            
            # Get labels and predictions
            y_test = true_labels[disease]
            valid_indices = y_test.isin([0, 1])
            
            if sum(valid_indices) == 0:
                print(f"Skipping {disease}: no valid labels")
                continue
                
            y_test = y_test[valid_indices]
            y_pred_proba = risk_predictions.loc[valid_indices, disease]
            
            if len(y_test) == 0 or len(set(y_test)) == 1:
                print(f"Skipping {disease}: insufficient class variety")
                continue
            
            # Calculate AUC
            try:
                auc_score = roc_auc_score(y_test, y_pred_proba)
            except:
                print(f"Failed to calculate AUC for {disease}")
                continue
            
            # Calculate ROC curve data
            fpr, tpr, roc_thresholds = roc_curve(y_test, y_pred_proba)
            
            # Calculate PR curve data
            precision, recall, pr_thresholds = precision_recall_curve(y_test, y_pred_proba)
            auprc = auc(recall, precision)
            
            # Find optimal threshold
            optimal_threshold, max_f1 = self.find_optimal_threshold(y_test, y_pred_proba)
            
            # Calculate binary predictions
            y_pred_binary = (y_pred_proba >= optimal_threshold).astype(int)
            
            # Calculate confusion matrix
            cm = confusion_matrix(y_test, y_pred_binary)
            tn, fp, fn, tp = cm.ravel() if cm.size == 4 else (0, 0, 0, 0)
            
            # Calculate additional metrics
            sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
            specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
            precision_metric = tp / (tp + fp) if (tp + fp) > 0 else 0
            
            # Calculate F1 scores
            f1_micro = f1_score(y_test, y_pred_binary, average='micro', zero_division=0)
            f1_macro = f1_score(y_test, y_pred_binary, average='macro', zero_division=0)
            f1_weighted = f1_score(y_test, y_pred_binary, average='weighted', zero_division=0)
            
            # Store AUC results
            auc_results.append({
                "Disease": disease,
                "AUC": auc_score,
                "AUPRC": auprc,
                "Sensitivity": sensitivity,
                "Specificity": specificity,
                "Precision": precision_metric,
                "F1_optimal": max_f1,
                "F1_micro": f1_micro,
                "F1_macro": f1_macro,
                "F1_weighted": f1_weighted,
                "Optimal_threshold": optimal_threshold,
                "Positive_samples": sum(y_test),
                "Total_samples": len(y_test)
            })
            
            # Store ROC curve data
            for fp_val, tp_val, thresh in zip(fpr, tpr, roc_thresholds):
                roc_curve_data.append({
                    "Disease": disease,
                    "FPR": fp_val,
                    "TPR": tp_val,
                    "Threshold": thresh,
                    "Type": "ROC"
                })
            
            # Store PR curve data
            for i, (prec, rec, thresh) in enumerate(zip(precision[:-1], recall[:-1], pr_thresholds)):
                pr_curve_data.append({
                    "Disease": disease,
                    "Precision": prec,
                    "Recall": rec,
                    "Threshold": thresh,
                    "Type": "PR"
                })
            
            # Add last point of PR curve
            pr_curve_data.append({
                "Disease": disease,
                "Precision": precision[-1],
                "Recall": recall[-1],
                "Threshold": np.nan,
                "Type": "PR"
            })
            
            print(f"Processed {disease}: AUC={auc_score:.4f}, F1={max_f1:.4f}")
        
        # Create DataFrames
        auc_df = pd.DataFrame(auc_results)
        roc_df = pd.DataFrame(roc_curve_data)
        pr_df = pd.DataFrame(pr_curve_data)
        
        # Save results if output directory is provided
        if output_dir:
            auc_df.to_csv(f"{output_dir}/auc_summary.csv", index=False)
            roc_df.to_csv(f"{output_dir}/roc_curve_data.csv", index=False)
            pr_df.to_csv(f"{output_dir}/pr_curve_data.csv", index=False)
            print(f"Results saved to {output_dir}")
        
        return auc_df, roc_df, pr_df