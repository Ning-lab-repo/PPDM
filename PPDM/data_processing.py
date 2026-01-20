"""
Data processing module for disease risk prediction
Handles loading, cleaning, and preprocessing of multiple data sources
"""

import pandas as pd
import numpy as np
from sklearn.impute import KNNImputer

class DataProcessor:
    def __init__(self):
        """Initialize DataProcessor with default settings"""
        self.train_country_code = 1  # England
        self.test_country_codes = [2, 3]  # Scotland and Wales
        
    def load_data(self, data_dir='../data'):
        """Load all required data files from the data directory"""
        data_paths = {
            'diagnoses': f'{data_dir}/all_icd_f_group_1000_10y_wear.csv',
            'movement': f'{data_dir}/ml/all_movement_features.csv',
            'validation': f'{data_dir}/participant_country_classification.csv',
            'baseline': f'{data_dir}/Baseline_characteristics_6.csv',
            'static': f'{data_dir}/ml/all_static_features.csv',
            'chronic': f'{data_dir}/chronic_diseases_before.csv'
        }
        
        data_dict = {}
        for name, path in data_paths.items():
            data_dict[name] = pd.read_csv(path)
            print(f"Loaded {name}: {data_dict[name].shape}")
        
        return data_dict
    
    def preprocess_data(self, data_dict):
        """Preprocess and merge all data sources"""
        # Extract validation info
        validation = data_dict['validation'][['Participant ID', 'Country']]
        country_mapping = {'England': 1, 'Scotland': 2, 'Wales': 3}
        validation['test'] = validation['Country'].map(country_mapping)
        
        # Get valid IDs from movement features
        valid_ids = data_dict['movement']['Participant ID'].unique()
        
        # Filter diagnoses data
        diagnoses = data_dict['diagnoses'].copy()
        diagnoses = diagnoses[diagnoses['Participant ID'].isin(valid_ids)]
        
        # Merge datasets with validation info
        datasets = {}
        
        # Diagnoses data
        datasets['diagnoses'] = pd.merge(
            validation[["Participant ID", 'test']], 
            diagnoses, 
            on='Participant ID', 
            how='right'
        )
        
        # Movement features
        datasets['movement'] = pd.merge(
            validation[["Participant ID", 'test']], 
            data_dict['movement'], 
            on='Participant ID', 
            how='right'
        )
        
        # Static features (with chronic diseases)
        static_features = pd.merge(
            data_dict['chronic'], 
            data_dict['static'], 
            on='Participant ID', 
            how='right'
        )
        datasets['static'] = pd.merge(
            validation[["Participant ID", 'test']], 
            static_features, 
            on='Participant ID', 
            how='right'
        )
        
        # Sort all datasets by Participant ID
        for name in datasets:
            datasets[name] = datasets[name].sort_values(by='Participant ID').reset_index(drop=True)
        
        return datasets
    
    def impute_missing_values(self, datasets):
        """Apply KNN imputation to handle missing values"""
        knn_imputer = KNNImputer(n_neighbors=5)
        
        # Impute movement features (from column 2 onward)
        movement_to_impute = datasets['movement'].iloc[:, 2:]
        datasets['movement'].iloc[:, 2:] = knn_imputer.fit_transform(movement_to_impute)
        
        # Impute static features (from column 2 onward)
        static_to_impute = datasets['static'].iloc[:, 2:]
        datasets['static'].iloc[:, 2:] = knn_imputer.fit_transform(static_to_impute)
        
        return datasets
    
    def split_data(self, datasets):
        """Split data into training and testing sets"""
        # Time gap feature
        time_gap = datasets['static'][["Participant ID", 'test', 'Time difference_x']].copy()
        datasets['static'] = datasets['static'].drop(columns=['Time difference_x'])
        
        # Split diagnoses
        train_diagnoses = datasets['diagnoses'][datasets['diagnoses']['test'] == self.train_country_code]
        test_diagnoses = datasets['diagnoses'][datasets['diagnoses']['test'].isin(self.test_country_codes)]
        
        # Split movement features and reshape
        train_movement = datasets['movement'].iloc[:, 2:][datasets['movement']['test'] == self.train_country_code]
        test_movement = datasets['movement'].iloc[:, 2:][datasets['movement']['test'].isin(self.test_country_codes)]
        
        train_movement = train_movement.values.reshape(-1, 48, 4)
        test_movement = test_movement.values.reshape(-1, 48, 4)
        
        # Split static features
        train_static = datasets['static'].iloc[:, 2:][datasets['static']['test'] == self.train_country_code]
        test_static = datasets['static'].iloc[:, 2:][datasets['static']['test'].isin(self.test_country_codes)]
        
        # Split time gap
        train_time_gap = time_gap.iloc[:, 2:][time_gap['test'] == self.train_country_code]
        test_time_gap = time_gap.iloc[:, 2:][time_gap['test'].isin(self.test_country_codes)]
        
        # Convert to numpy arrays
        train_static = train_static.to_numpy()
        test_static = test_static.to_numpy()
        train_time_gap = train_time_gap.to_numpy()
        test_time_gap = test_time_gap.to_numpy()
        
        # Get test participant IDs
        test_participant_ids = datasets['movement']['Participant ID'][
            datasets['movement']['test'].isin(self.test_country_codes)
        ]
        
        split_data = {
            'train_diagnoses': train_diagnoses,
            'test_diagnoses': test_diagnoses,
            'train_movement': train_movement,
            'test_movement': test_movement,
            'train_static': train_static,
            'test_static': test_static,
            'train_time_gap': train_time_gap,
            'test_time_gap': test_time_gap,
            'test_participant_ids': test_participant_ids
        }
        
        return split_data