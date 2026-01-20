"""
Main execution script for optimal pattern generation
Coordinates data loading, pattern generation, and output
"""

import pandas as pd
import numpy as np
import argparse
from datetime import datetime

# Import custom modules
from config import PATHS, ANALYSIS, GROUPS
from data_loader import DataLoader
from pattern_generator import PatternGenerator
from utils import OutputUtils

def main(args):
    """Main execution pipeline"""
    print("="*60)
    print("OPTIMAL MOVEMENT PATTERN GENERATOR")
    print("="*60)
    
    # Initialize modules
    data_loader = DataLoader()
    pattern_gen = PatternGenerator()
    output_utils = OutputUtils()
    
    # Create output directory
    output_utils.create_output_directory()
    
    # Step 1: Load and preprocess data
    print("\n1. LOADING AND PREPROCESSING DATA")
    print("-"*40)
    
    score_data, feature_data, number_data, _, id_col = data_loader.load_and_align_data()
    score_filtered, number_filtered, icd_list = data_loader.filter_icd_codes(score_data, number_data)
    reshaped_data, feature_columns = data_loader.reshape_movement_features(feature_data)
    
    # Step 2: Generate optimal patterns for each group
    print("\n2. GENERATING OPTIMAL PATTERNS")
    print("-"*40)
    
    all_patterns = {}
    
    # Generate patterns for each stratification group
    for group in GROUPS:
        group_name = group.replace('_', ' ').title()
        if group.startswith('age_under'):
            group_name = f"Age Under {ANALYSIS['age_threshold']}"
        elif group.startswith('age_') and '_plus' in group:
            group_name = f"Age {ANALYSIS['age_threshold']} Plus"
        
        print(f"\nProcessing group: {group_name}")
        
        # Get stratified data
        stratified_data = data_loader.get_stratified_data(feature_data, group)
        print(f"  Sample size: {len(stratified_data)} participants")
        
        # Generate optimal pattern for this group
        group_pattern = pattern_gen.generate_group_optimal_pattern(
            score_data=score_filtered,
            feature_data=stratified_data,
            reshaped_data=reshaped_data,
            icd_list=icd_list,
            group_type=group_name
        )
        
        all_patterns[group_name] = group_pattern
    
    # Step 3: Generate disease-specific patterns (optional)
    print("\n3. GENERATING DISEASE-SPECIFIC PATTERNS")
    print("-"*40)
    
    disease_patterns = {}
    if args.generate_disease_specific:
        print(f"Generating patterns for {min(args.max_diseases, len(icd_list))} diseases...")
        
        for i, icd in enumerate(icd_list[:args.max_diseases]):
            if i % 50 == 0:
                print(f"  Processed {i} diseases...")
            
            # Generate pattern for this disease
            risk_values = score_filtered[icd].values
            threshold_low = np.percentile(risk_values, 1)
            threshold_high = np.percentile(risk_values, 99)
            
            low_risk_idx = np.where(risk_values <= threshold_low)[0]
            high_risk_idx = np.where(risk_values >= threshold_high)[0]
            
            if len(low_risk_idx) >= ANALYSIS['min_sample_size'] and \
               len(high_risk_idx) >= ANALYSIS['min_sample_size']:
                
                low_risk_pattern = reshaped_data[low_risk_idx].mean(axis=0)
                high_risk_pattern = reshaped_data[high_risk_idx].mean(axis=0)
                
                # Calculate optimal direction
                risk_difference = high_risk_pattern - low_risk_pattern
                optimal_direction = -risk_difference
                
                # Limit adjustment
                adjustment_magnitude = np.abs(optimal_direction).max()
                if adjustment_magnitude > 0:
                    scaling_factor = min(ANALYSIS['max_adjustment'] / adjustment_magnitude, 1.0)
                    optimal_direction = optimal_direction * scaling_factor
                
                disease_pattern = low_risk_pattern + optimal_direction * 1
                disease_pattern = np.maximum(disease_pattern, 0)
                
                # Normalize
                row_sums = disease_pattern.sum(axis=1, keepdims=True)
                disease_pattern = disease_pattern / row_sums
                
                # Apply constraints
                baseline_pattern = reshaped_data.mean(axis=0)
                disease_pattern = pattern_gen.apply_temporal_constraints(disease_pattern, baseline_pattern)
                
                disease_patterns[icd] = disease_pattern
        
        print(f"Generated {len(disease_patterns)} disease-specific patterns")
    
    # Step 4: Save all patterns
    print("\n4. SAVING RESULTS")
    print("-"*40)
    
    saved_files = []
    
    # Save group patterns
    for pattern_name, pattern_data in all_patterns.items():
        output_path = output_utils.save_pattern_as_csv(
            pattern_data, pattern_name, feature_columns
        )
        saved_files.append(output_path)
        print(f"✓ Saved: {pattern_name} pattern")
    
    # Save disease-specific patterns if generated
    if disease_patterns:
        disease_patterns_list = []
        for icd, pattern_data in disease_patterns.items():
            flattened_pattern = pattern_data.reshape(-1)
            disease_patterns_list.append(flattened_pattern)
        
        disease_patterns_df = pd.DataFrame(disease_patterns_list, columns=feature_columns)
        disease_patterns_df['Disease_ICD'] = list(disease_patterns.keys())
        disease_patterns_df['Pattern_Name'] = disease_patterns_df['Disease_ICD'].apply(
            lambda x: f'Optimal_{x}'
        )
        
        disease_output_path = os.path.join(
            PATHS['output_dir'], 'disease_specific_patterns.csv'
        )
        disease_patterns_df.to_csv(disease_output_path, index=False)
        saved_files.append(disease_output_path)
        print(f"✓ Saved: Disease-specific patterns ({len(disease_patterns)} diseases)")
    
    # Step 5: Generate summary statistics
    print("\n5. SUMMARY STATISTICS")
    print("-"*40)
    
    output_utils.generate_summary_statistics(all_patterns)
    output_utils.print_execution_summary(all_patterns, len(disease_patterns))
    
    print("\n" + "="*60)
    print("EXECUTION COMPLETED SUCCESSFULLY")
    print("="*60)

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Optimal Movement Pattern Generator")
    
    parser.add_argument("--generate_disease_specific", action="store_true",
                       help="Generate disease-specific patterns (can be time-consuming)")
    parser.add_argument("--max_diseases", type=int, default=100,
                       help="Maximum number of disease-specific patterns to generate")
    parser.add_argument("--data_dir", type=str, default="./data",
                       help="Directory containing input data")
    
    args = parser.parse_args()
    
    # Update paths based on command line argument
    if args.data_dir != "./data":
        PATHS['score_file'] = f"{args.data_dir}/predictions/risk_predictions.csv"
        PATHS['feature_file'] = f"{args.data_dir}/features/movement_features.csv"
        PATHS['number_file'] = f"{args.data_dir}/diagnoses/disease_counts.csv"
        PATHS['age_sex_file'] = f"{args.data_dir}/demographics/age_sex.csv"
    
    # Run main pipeline
    main(args)