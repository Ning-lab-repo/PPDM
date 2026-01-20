"""
Configuration file for optimal pattern generation
All paths and parameters are centralized here
"""

# ======================================
# File Paths (Virtual paths for GitHub)
# ======================================

PATHS = {
    'score_file': './data/predictions/risk_predictions.csv',
    'feature_file': './data/features/movement_features.csv',
    'number_file': './data/diagnoses/disease_counts.csv',
    'age_sex_file': './data/demographics/age_sex.csv',
    'output_dir': './results/optimal_patterns'
}

# ======================================
# Analysis Parameters
# ======================================

ANALYSIS = {
    'icd_filter': 'N',  # Keep ICD codes with first letter <= 'N'
    'age_threshold': 60,  # Age threshold for stratification
    'risk_percentile_low': 5,  # Percentile for low risk group
    'risk_percentile_high': 95,  # Percentile for high risk group
    'min_sample_size': 10,  # Minimum sample size for pattern generation
    'max_adjustment': 1.2,  # Maximum adjustment factor for patterns
    'smoothing_factor': 0.0  # Smoothing between optimal and baseline
}

# ======================================
# Temporal Constraints
# ======================================

CONSTRAINTS = {
    'night_hours': list(range(0, 6)) + list(range(22, 24)),  # 22:00-6:00
    'day_hours': list(range(8, 20)),  # 8:00-20:00
    'max_night_activity': 0.1,  # Maximum activity during night
    'max_day_sedentary': 0.8,  # Maximum sedentary during day
}

# ======================================
# Activity States
# ======================================

STATES = {
    'names': ['Sedentary', 'Light', 'Moderate-Vigorous', 'Sleep'],
    'indices': [0, 1, 2, 3]
}

# ======================================
# Stratification Groups
# ======================================

GROUPS = [
    'all',
    'male',
    'female',
    f'age_under{ANALYSIS["age_threshold"]}',
    f'age_{ANALYSIS["age_threshold"]}_plus'
]