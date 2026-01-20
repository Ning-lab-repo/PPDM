"""
Optimal pattern generation algorithms
Implements various strategies for finding optimal movement patterns
"""

import numpy as np
import pandas as pd
from config import ANALYSIS, CONSTRAINTS, STATES

class PatternGenerator:
    def __init__(self, params=None):
        """Initialize PatternGenerator with parameters"""
        self.params = params if params else ANALYSIS
        self.constraints = CONSTRAINTS
        self.states = STATES
    
    def get_low_risk_patterns(self, score_data, reshaped_data, icd_list, stratified_ids=None):
        """Get patterns from low-risk individuals"""
        direct_patterns = {}
        
        if stratified_ids is not None:
            score_data = score_data[score_data[score_data.columns[0]].isin(stratified_ids)]
            reshaped_indices = [i for i, pid in enumerate(score_data[score_data.columns[0]]) 
                               if pid in stratified_ids]
            reshaped_data = reshaped_data[reshaped_indices]
        
        for icd in icd_list:
            risk_values = score_data[icd].values
            threshold = np.percentile(risk_values, self.params['risk_percentile_low'])
            low_risk_idx = np.where(risk_values <= threshold)[0]
            
            if len(low_risk_idx) >= self.params['min_sample_size']:
                low_risk_pattern = reshaped_data[low_risk_idx].mean(axis=0)
                # Normalize to ensure probabilities sum to 1
                row_sums = low_risk_pattern.sum(axis=1, keepdims=True)
                low_risk_pattern = low_risk_pattern / row_sums
                direct_patterns[icd] = low_risk_pattern
        
        return direct_patterns
    
    def get_contrastive_patterns(self, score_data, reshaped_data, icd_list, stratified_ids=None):
        """Get patterns using contrastive learning approach"""
        contrastive_patterns = {}
        
        if stratified_ids is not None:
            score_data = score_data[score_data[score_data.columns[0]].isin(stratified_ids)]
            reshaped_indices = [i for i, pid in enumerate(score_data[score_data.columns[0]]) 
                               if pid in stratified_ids]
            reshaped_data = reshaped_data[reshaped_indices]
        
        for icd in icd_list:
            risk_values = score_data[icd].values
            threshold_low = np.percentile(risk_values, 1)  # Bottom 1%
            threshold_high = np.percentile(risk_values, 99)  # Top 1%
            
            low_risk_idx = np.where(risk_values <= threshold_low)[0]
            high_risk_idx = np.where(risk_values >= threshold_high)[0]
            
            if (len(low_risk_idx) < self.params['min_sample_size'] or 
                len(high_risk_idx) < self.params['min_sample_size']):
                continue
            
            low_risk_pattern = reshaped_data[low_risk_idx].mean(axis=0)
            high_risk_pattern = reshaped_data[high_risk_idx].mean(axis=0)
            
            # Calculate optimal direction (away from high risk)
            risk_difference = high_risk_pattern - low_risk_pattern
            optimal_direction = -risk_difference
            
            # Limit adjustment magnitude
            adjustment_magnitude = np.abs(optimal_direction).max()
            if adjustment_magnitude > 0:
                scaling_factor = min(self.params['max_adjustment'] / adjustment_magnitude, 1.0)
                optimal_direction = optimal_direction * scaling_factor
            
            # Apply adjustment
            adjusted_pattern = low_risk_pattern + optimal_direction * 1
            adjusted_pattern = np.maximum(adjusted_pattern, 0)
            
            # Normalize
            row_sums = adjusted_pattern.sum(axis=1, keepdims=True)
            adjusted_pattern = adjusted_pattern / row_sums
            
            contrastive_patterns[icd] = adjusted_pattern
        
        return contrastive_patterns
    
    def apply_temporal_constraints(self, pattern, baseline_pattern=None):
        """Apply temporal feasibility constraints to pattern"""
        constrained_pattern = pattern.copy()
        
        # Night time constraints
        for hour in self.constraints['night_hours']:
            if constrained_pattern[hour, self.states['indices'][2]] > self.constraints['max_night_activity']:
                reduction = constrained_pattern[hour, self.states['indices'][2]] - self.constraints['max_night_activity']
                constrained_pattern[hour, self.states['indices'][2]] = self.constraints['max_night_activity']
                # Redistribute to sleep and light activity
                constrained_pattern[hour, self.states['indices'][3]] += reduction * 0.7
                constrained_pattern[hour, self.states['indices'][1]] += reduction * 0.3
        
        # Day time constraints
        for hour in self.constraints['day_hours']:
            if constrained_pattern[hour, self.states['indices'][0]] > self.constraints['max_day_sedentary']:
                reduction = constrained_pattern[hour, self.states['indices'][0]] - self.constraints['max_day_sedentary']
                constrained_pattern[hour, self.states['indices'][0]] = self.constraints['max_day_sedentary']
                # Redistribute to light and moderate-vigorous activity
                constrained_pattern[hour, self.states['indices'][1]] += reduction * 0.6
                constrained_pattern[hour, self.states['indices'][2]] += reduction * 0.4
        
        # Ensure rows sum to 1
        row_sums = constrained_pattern.sum(axis=1, keepdims=True)
        constrained_pattern = constrained_pattern / row_sums
        
        return constrained_pattern
    
    def generate_group_optimal_pattern(self, score_data, feature_data, reshaped_data, 
                                       icd_list, group_type="Overall"):
        """Generate optimal pattern for a specific group"""
        print(f"\nGenerating {group_type} optimal pattern...")
        
        # Get stratified data if not overall
        if group_type != "Overall":
            stratified_data = feature_data
            stratified_ids = set(stratified_data['Participant ID'])
        else:
            stratified_data = feature_data
            stratified_ids = None
        
        # Generate patterns using different methods
        direct_patterns = self.get_low_risk_patterns(
            score_data, reshaped_data, icd_list, stratified_ids
        )
        contrastive_patterns = self.get_contrastive_patterns(
            score_data, reshaped_data, icd_list, stratified_ids
        )
        
        print(f"  Direct patterns: {len(direct_patterns)} diseases")
        print(f"  Contrastive patterns: {len(contrastive_patterns)} diseases")
        
        # Combine available patterns (prefer contrastive)
        available_patterns = {}
        for icd in icd_list:
            if icd in contrastive_patterns:
                available_patterns[icd] = contrastive_patterns[icd]
            elif icd in direct_patterns:
                available_patterns[icd] = direct_patterns[icd]
        
        print(f"  Available patterns: {len(available_patterns)} diseases")
        
        # Weighted average based on disease prevalence
        if len(available_patterns) == 0:
            print(f"  Warning: No patterns available, using global low-risk pattern")
            if stratified_ids is not None:
                stratified_score = score_data[score_data[score_data.columns[0]].isin(stratified_ids)]
                threshold_global = np.percentile(stratified_score[icd_list].mean(axis=1), 10)
                global_low_risk = np.where(stratified_score[icd_list].mean(axis=1) <= threshold_global)[0]
                stratified_reshaped = reshaped_data[[i for i, pid in enumerate(
                    score_data[score_data.columns[0]]) if pid in stratified_ids]]
                weighted_optimal = stratified_reshaped[global_low_risk].mean(axis=0)
            else:
                # Fallback to overall baseline
                weighted_optimal = reshaped_data.mean(axis=0)
        else:
            # Calculate weights based on disease prevalence
            case_counts = {icd: len(available_patterns[icd]) for icd in available_patterns}
            total_cases = sum(case_counts.values())
            weights = {icd: count/total_cases for icd, count in case_counts.items()}
            
            # Compute weighted average
            weighted_optimal = np.zeros((48, 4))
            for icd, pattern in available_patterns.items():
                weighted_optimal += pattern * weights[icd]
            
            # Normalize
            row_sums = weighted_optimal.sum(axis=1, keepdims=True)
            weighted_optimal = weighted_optimal / row_sums
        
        # Apply temporal constraints
        if stratified_ids is not None:
            stratified_reshaped = reshaped_data[[i for i, pid in enumerate(
                score_data[score_data.columns[0]]) if pid in stratified_ids]]
            baseline_pattern = stratified_reshaped.mean(axis=0)
        else:
            baseline_pattern = reshaped_data.mean(axis=0)
        
        final_optimal = self.apply_temporal_constraints(weighted_optimal, baseline_pattern)
        
        # Smooth with baseline if needed
        if self.params['smoothing_factor'] > 0:
            final_optimal = (1 - self.params['smoothing_factor']) * final_optimal + \
                           self.params['smoothing_factor'] * baseline_pattern
        
        # Final normalization
        row_sums = final_optimal.sum(axis=1, keepdims=True)
        final_optimal = final_optimal / row_sums
        
        return final_optimal