"""
Utility functions for data processing and output
"""

import os
import pandas as pd
import numpy as np
from config import PATHS, STATES

class OutputUtils:
    def __init__(self):
        """Initialize OutputUtils"""
        self.paths = PATHS
        self.states = STATES
    
    def create_output_directory(self):
        """Create output directory if it doesn't exist"""
        os.makedirs(self.paths['output_dir'], exist_ok=True)
        print(f"Output directory created: {self.paths['output_dir']}")
    
    def save_pattern_as_csv(self, pattern_data, pattern_name, feature_columns):
        """Save pattern as CSV file"""
        # Flatten pattern to 1D array
        flattened_pattern = pattern_data.reshape(-1)
        
        # Create DataFrame
        pattern_df = pd.DataFrame([flattened_pattern], columns=feature_columns)
        pattern_df['Pattern_Name'] = pattern_name
        
        # Save to file
        output_path = os.path.join(
            self.paths['output_dir'], 
            f"{pattern_name.lower().replace(' ', '_')}_movement_pattern.csv"
        )
        pattern_df.to_csv(output_path, index=False)
        
        return output_path
    
    def generate_summary_statistics(self, patterns_dict):
        """Generate summary statistics for patterns"""
        print("\n" + "="*60)
        print("PATTERN SUMMARY STATISTICS")
        print("="*60)
        
        for pattern_name, pattern_data in patterns_dict.items():
            print(f"\n{pattern_name}:")
            for i, state_name in enumerate(self.states['names']):
                percentage = pattern_data[:, i].mean() * 100
                print(f"  {state_name}: {percentage:.1f}%")
        
        return self._compare_age_groups(patterns_dict)
    
    def _compare_age_groups(self, patterns_dict):
        """Compare age group patterns"""
        age_key = f'age_under{60}'  # Using 60 as threshold from config
        age_plus_key = f'age_{60}_plus'
        
        if age_key in patterns_dict and age_plus_key in patterns_dict:
            print(f"\nAGE GROUP COMPARISON (threshold: {60} years):")
            print("-"*40)
            
            under_data = patterns_dict[age_key]
            plus_data = patterns_dict[age_plus_key]
            
            for i, state_name in enumerate(self.states['names']):
                under_pct = under_data[:, i].mean() * 100
                plus_pct = plus_data[:, i].mean() * 100
                diff = plus_pct - under_pct
                print(f"  {state_name}: Under {60} = {under_pct:.1f}% | "
                      f"{60}+ = {plus_pct:.1f}% | Difference = {diff:+.1f}%")
    
    def print_execution_summary(self, patterns_dict, disease_patterns_count):
        """Print final execution summary"""
        print("\n" + "="*60)
        print("EXECUTION SUMMARY")
        print("="*60)
        print(f"✓ Generated overall patterns: {len(patterns_dict)}")
        print(f"✓ Generated disease-specific patterns: {disease_patterns_count}")
        print(f"✓ All results saved to: {self.paths['output_dir']}")
        print("="*60)