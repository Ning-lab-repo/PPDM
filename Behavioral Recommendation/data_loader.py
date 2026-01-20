"""
Data loading and preprocessing module
Handles data loading, filtering, and alignment
"""

import pandas as pd
import numpy as np
from config import PATHS, ANALYSIS

class DataLoader:
    def __init__(self):
        """Initialize DataLoader with configuration"""
        self.paths = PATHS
        self.params = ANALYSIS
        
    def load_and_align_data(self):
        """Load and align all required data sources"""
        print("Loading data from files...")
        
        # Load all data files
        score_data = pd.read_csv(self.paths['score_file'])
        feature_data = pd.read_csv(self.paths['feature_file'])
        number_data = pd.read_csv(self.paths['number_file'])
        age_sex_data = pd.read_csv(self.paths['age_sex_file'])
        
        print(f"Loaded score data: {score_data.shape}")
        print(f"Loaded feature data: {feature_data.shape}")
        print(f"Loaded number data: {number_data.shape}")
        print(f"Loaded demographics data: {age_sex_data.shape}")
        
        # Identify participant ID column
        id_col = score_data.columns[0]
        
        # Find common participants across all datasets
        common_ids = set(score_data[id_col]) & set(feature_data['Participant ID']) & \
                     set(number_data['Participant ID']) & set(age_sex_data['Participant ID'])
        
        print(f"Found {len(common_ids)} common participants")
        
        # Filter data to common participants
        score_data = score_data[score_data[id_col].isin(common_ids)].reset_index(drop=True)
        feature_data = feature_data[feature_data['Participant ID'].isin(common_ids)].reset_index(drop=True)
        number_data = number_data[number_data['Participant ID'].isin(common_ids)].reset_index(drop=True)
        age_sex_data = age_sex_data[age_sex_data['Participant ID'].isin(common_ids)].reset_index(drop=True)
        
        # Merge demographics with features
        feature_data = feature_data.merge(
            age_sex_data[['Participant ID', 'age', 'Sex']], 
            on='Participant ID'
        )
        
        return score_data, feature_data, number_data, age_sex_data, id_col
    
    def filter_icd_codes(self, score_data, number_data):
        """Filter ICD codes based on first letter criterion"""
        all_icds = [col for col in score_data.columns if col != score_data.columns[0]]
        filtered_icds = [
            icd for icd in all_icds 
            if isinstance(icd, str) and icd[0].upper() <= self.params['icd_filter']
        ]
        
        print(f"Filtered to {len(filtered_icds)} ICD codes (first letter <= {self.params['icd_filter']})")
        
        # Filter score and number data
        score_filtered = score_data[[score_data.columns[0]] + filtered_icds]
        number_filtered = number_data[['Participant ID'] + filtered_icds]
        
        return score_filtered, number_filtered, filtered_icds
    
    def reshape_movement_features(self, feature_data):
        """Reshape movement features to 3D array"""
        # Extract movement features (excluding ID and demographic columns)
        movement_columns = [col for col in feature_data.columns 
                          if col not in ['Participant ID', 'age', 'Sex']]
        movement_data = feature_data[movement_columns].values
        
        # Reshape to 3D: (participants, 48 hours, 4 states)
        n_participants = movement_data.shape[0]
        reshaped_data = movement_data.reshape(n_participants, 48, 4)
        
        return reshaped_data, movement_columns
    
    def get_stratified_data(self, feature_data, group_type='all'):
        """Get data stratified by demographic groups"""
        if group_type == 'all':
            return feature_data
        elif group_type == 'male':
            return feature_data[feature_data['Sex'] == 1]
        elif group_type == 'female':
            return feature_data[feature_data['Sex'] == 0]
        elif group_type == f'age_under{self.params["age_threshold"]}':
            return feature_data[feature_data['age'] < self.params['age_threshold']]
        elif group_type == f'age_{self.params["age_threshold"]}_plus':
            return feature_data[feature_data['age'] >= self.params['age_threshold']]
        else:
            raise ValueError(f"Unknown group type: {group_type}")