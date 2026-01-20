"""
Training pipeline for disease risk prediction models
Handles model training, validation, and parallel processing
"""

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, confusion_matrix
from joblib import Parallel, delayed

class ModelTrainer:
    def __init__(self, device='cuda' if torch.cuda.is_available() else 'cpu'):
        """Initialize ModelTrainer with device settings"""
        self.device = torch.device(device)
        print(f"Using device: {self.device}")
        
    def train_single_model(self, disease, train_data, model_params):
        """Train a single model for a specific disease"""
        # Unpack training data
        y = train_data['diagnoses'][disease]
        valid_indices = y.isin([0, 1])
        y = y[valid_indices]
        
        if len(y) < 10:
            print(f"Skipping {disease}: insufficient samples ({len(y)})")
            return None
            
        # Prepare features
        X = train_data['movement'][valid_indices]
        workday_data = X[:, :24, :]
        restday_data = X[:, 24:, :]
        
        # Convert to tensors
        workday_tensor = torch.tensor(workday_data, dtype=torch.float32)
        restday_tensor = torch.tensor(restday_data, dtype=torch.float32)
        static_tensor = torch.tensor(train_data['static'][valid_indices], dtype=torch.float32)
        time_gap_tensor = torch.tensor(train_data['time_gap'][valid_indices], dtype=torch.float32)
        y_tensor = torch.tensor(y.values, dtype=torch.float32).unsqueeze(1)
        
        # Split into train/validation sets
        workday_train, workday_val, restday_train, restday_val, static_train, static_val, \
        time_gap_train, time_gap_val, y_train, y_val = train_test_split(
            workday_tensor, restday_tensor, static_tensor, time_gap_tensor, y_tensor,
            test_size=0.2, random_state=42, stratify=y_tensor.numpy()
        )
        
        # Create data loaders
        train_dataset = TensorDataset(workday_train, restday_train, static_train, 
                                      time_gap_train, y_train)
        train_loader = DataLoader(train_dataset, batch_size=128, shuffle=True)
        
        # Initialize model
        from model_architecture import LSTM_SelfAttention_Model
        model = LSTM_SelfAttention_Model(
            input_size=4,
            hidden_size=32,
            output_size=1,
            static_feature_size=train_data['static'].shape[1],
            lstm_dropout_rate=0.5,
            attention_dropout_rate=0.3,
            use_time_aware=model_params.get('use_time_aware', True),
            use_dnn=model_params.get('use_dnn', True)
        ).to(self.device)
        
        # Training configuration
        criterion = nn.BCEWithLogitsLoss()
        optimizer = optim.Adam(model.parameters(), lr=0.001)
        
        # Training loop
        model.train()
        for epoch in range(20):
            epoch_loss = 0
            for batch in train_loader:
                workday_inputs, restday_inputs, static_inputs, time_gap_inputs, labels = batch
                
                # Move to device
                workday_inputs = workday_inputs.to(self.device)
                restday_inputs = restday_inputs.to(self.device)
                static_inputs = static_inputs.to(self.device)
                time_gap_inputs = time_gap_inputs.to(self.device)
                labels = labels.to(self.device)
                
                optimizer.zero_grad()
                outputs, _, _ = model(workday_inputs, restday_inputs, static_inputs, time_gap_inputs)
                loss = criterion(outputs, labels)
                loss.backward()
                optimizer.step()
                epoch_loss += loss.item()
                
            print(f"  Disease: {disease}, Epoch {epoch+1}/20, Loss: {epoch_loss/len(train_loader):.4f}")
        
        # Evaluate on validation set
        model.eval()
        with torch.no_grad():
            workday_val = workday_val.to(self.device)
            restday_val = restday_val.to(self.device)
            static_val = static_val.to(self.device)
            time_gap_val = time_gap_val.to(self.device)
            
            y_pred_proba, workday_attn_weights, restday_attn_weights = model(
                workday_val, restday_val, static_val, time_gap_val
            )
            y_pred_proba = y_pred_proba.squeeze()
            y_pred = (torch.sigmoid(y_pred_proba) > 0.5).float()
            
            auc = roc_auc_score(y_val.cpu().numpy(), torch.sigmoid(y_pred_proba).cpu().numpy())
            cm = confusion_matrix(y_val.cpu().numpy(), y_pred.cpu().numpy())
            TP = cm[1, 1]
            FN = cm[1, 0]
            sensitivity = TP / (TP + FN) if (TP + FN) != 0 else 0
        
        print(f"  Disease: {disease}, AUC: {auc:.4f}, Sensitivity: {sensitivity:.4f}")
        
        return {
            'disease': disease,
            'model': model.cpu(),  # Move to CPU for saving
            'auc': auc,
            'sensitivity': sensitivity,
            'attention_weights': {
                'workday': workday_attn_weights,
                'restday': restday_attn_weights
            }
        }
    
    def train_all_models(self, train_data, model_params, n_jobs=1):
        """Train models for all diseases in parallel"""
        disease_columns = train_data['diagnoses'].columns[2:]  # Skip Participant ID and test columns
        
        print(f"Training models for {len(disease_columns)} diseases")
        print(f"Using {n_jobs} parallel jobs")
        
        # Prepare training data dict
        train_data_dict = {
            'diagnoses': train_data['diagnoses'],
            'movement': train_data['movement'],
            'static': train_data['static'],
            'time_gap': train_data['time_gap']
        }
        
        # Parallel training
        if n_jobs > 1:
            results = Parallel(n_jobs=n_jobs)(
                delayed(self.train_single_model)(disease, train_data_dict, model_params)
                for disease in disease_columns
            )
        else:
            # Sequential training for debugging
            results = []
            for disease in disease_columns:
                result = self.train_single_model(disease, train_data_dict, model_params)
                results.append(result)
        
        # Filter out None results (diseases with insufficient data)
        results = [r for r in results if r is not None]
        
        # Collect results
        trained_models = {}
        metrics = {}
        
        for result in results:
            disease = result['disease']
            trained_models[disease] = result['model']
            metrics[disease] = {
                'auc': result['auc'],
                'sensitivity': result['sensitivity']
            }
        
        print(f"\nTraining completed for {len(trained_models)} diseases")
        
        # Print summary
        print("\nSummary of results:")
        for disease, metric in metrics.items():
            print(f"  {disease}: AUC={metric['auc']:.4f}, Sensitivity={metric['sensitivity']:.4f}")
        
        return trained_models, metrics