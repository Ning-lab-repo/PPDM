"""
Main execution script for disease risk prediction pipeline
Orchestrates data processing, model training, and evaluation
"""

import argparse
import joblib
import os
from datetime import datetime

# Import custom modules
from data_processing import DataProcessor
from training_pipeline import ModelTrainer
from evaluation import ModelEvaluator

def main(args):
    """Main pipeline execution"""
    print("=" * 60)
    print("DISEASE RISK PREDICTION PIPELINE")
    print("=" * 60)
    
    # Create output directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = f"./results/{timestamp}"
    os.makedirs(output_dir, exist_ok=True)
    
    # Step 1: Data Processing
    print("\n1. LOADING AND PROCESSING DATA")
    print("-" * 40)
    
    processor = DataProcessor()
    
    # Load data
    data_dict = processor.load_data(data_dir=args.data_dir)
    
    # Preprocess data
    datasets = processor.preprocess_data(data_dict)
    
    # Impute missing values
    datasets = processor.impute_missing_values(datasets)
    
    # Split data
    split_data = processor.split_data(datasets)
    
    # Save processed data
    joblib.dump(split_data, f"{output_dir}/processed_data.pkl")
    print(f"Processed data saved to {output_dir}/processed_data.pkl")
    
    # Step 2: Model Training
    print("\n2. TRAINING MODELS")
    print("-" * 40)
    
    if args.train_models:
        trainer = ModelTrainer(device=args.device)
        
        # Prepare training data
        train_data = {
            'diagnoses': split_data['train_diagnoses'],
            'movement': split_data['train_movement'],
            'static': split_data['train_static'],
            'time_gap': split_data['train_time_gap']
        }
        
        # Model parameters
        model_params = {
            'use_time_aware': args.use_time_aware,
            'use_dnn': args.use_dnn
        }
        
        # Train all models
        trained_models, metrics = trainer.train_all_models(
            train_data, 
            model_params, 
            n_jobs=args.n_jobs
        )
        
        # Save models
        models_dir = f"{output_dir}/models"
        os.makedirs(models_dir, exist_ok=True)
        
        for disease, model in trained_models.items():
            model_path = f"{models_dir}/{disease}_model.joblib"
            joblib.dump(model, model_path)
            print(f"Saved model for {disease}")
        
        # Save metrics
        metrics_df = pd.DataFrame(metrics).T
        metrics_df.to_csv(f"{output_dir}/training_metrics.csv")
        print(f"\nTraining metrics saved to {output_dir}/training_metrics.csv")
    else:
        print("Skipping model training (--train_models flag not set)")
        models_dir = args.model_dir
    
    # Step 3: Model Evaluation
    print("\n3. EVALUATING MODELS")
    print("-" * 40)
    
    evaluator = ModelEvaluator(device=args.device)
    
    # Load models (either trained or pre-trained)
    if args.train_models:
        models = trained_models
    else:
        models = evaluator.load_models(models_dir)
    
    if not models:
        print("No models found for evaluation")
        return
    
    # Generate predictions
    print(f"Generating predictions for {len(models)} diseases...")
    
    risk_predictions = evaluator.predict_in_batches(
        models=models,
        movement_data=split_data['test_movement'],
        participant_ids=split_data['test_participant_ids'],
        batch_size=args.batch_size
    )
    
    # Save predictions
    predictions_path = f"{output_dir}/risk_predictions.csv"
    risk_predictions.to_csv(predictions_path, index=False)
    print(f"Predictions saved to {predictions_path}")
    
    # Comprehensive evaluation
    print("\nPerforming comprehensive evaluation...")
    auc_df, roc_df, pr_df = evaluator.evaluate_all_metrics(
        risk_predictions=risk_predictions,
        true_labels=split_data['test_diagnoses'],
        output_dir=f"{output_dir}/evaluation"
    )
    
    # Print summary statistics
    print("\n4. SUMMARY STATISTICS")
    print("-" * 40)
    print(f"Number of diseases evaluated: {len(auc_df)}")
    print(f"Average AUC: {auc_df['AUC'].mean():.4f}")
    print(f"Average AUPRC: {auc_df['AUPRC'].mean():.4f}")
    print(f"Average Sensitivity: {auc_df['Sensitivity'].mean():.4f}")
    print(f"Average Specificity: {auc_df['Specificity'].mean():.4f}")
    
    # Save final summary
    summary_path = f"{output_dir}/pipeline_summary.txt"
    with open(summary_path, 'w') as f:
        f.write("Disease Risk Prediction Pipeline Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Timestamp: {timestamp}\n")
        f.write(f"Data directory: {args.data_dir}\n")
        f.write(f"Models trained: {args.train_models}\n")
        f.write(f"Number of diseases: {len(auc_df)}\n\n")
        f.write("Average Metrics:\n")
        f.write(f"  AUC: {auc_df['AUC'].mean():.4f}\n")
        f.write(f"  AUPRC: {auc_df['AUPRC'].mean():.4f}\n")
        f.write(f"  Sensitivity: {auc_df['Sensitivity'].mean():.4f}\n")
        f.write(f"  Specificity: {auc_df['Specificity'].mean():.4f}\n")
    
    print(f"\nPipeline completed successfully!")
    print(f"All results saved to: {output_dir}")
    print("=" * 60)

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Disease Risk Prediction Pipeline")
    
    # Data arguments
    parser.add_argument("--data_dir", type=str, default="../data",
                       help="Directory containing input data files")
    
    # Model arguments
    parser.add_argument("--model_dir", type=str, default="./models",
                       help="Directory containing pre-trained models")
    parser.add_argument("--train_models", action="store_true",
                       help="Train new models (if not set, use pre-trained models)")
    parser.add_argument("--use_time_aware", action="store_true", default=True,
                       help="Use time-aware module in the model")
    parser.add_argument("--use_dnn", action="store_true", default=True,
                       help="Use DNN module for static features")
    
    # Training arguments
    parser.add_argument("--n_jobs", type=int, default=4,
                       help="Number of parallel jobs for training")
    parser.add_argument("--device", type=str, default="cuda",
                       help="Device for training (cuda/cpu)")
    
    # Evaluation arguments
    parser.add_argument("--batch_size", type=int, default=5000,
                       help="Batch size for prediction")
    
    args = parser.parse_args()
    
    # Run main pipeline
    main(args)