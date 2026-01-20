"""
Neural network architecture definition
LSTM with self-attention for disease risk prediction
"""

import torch
import torch.nn as nn

class LSTM_SelfAttention_Model(nn.Module):
    """LSTM model with self-attention mechanism for time series classification"""
    
    def __init__(self, input_size, hidden_size, output_size, static_feature_size, 
                 lstm_dropout_rate=0.5, attention_dropout_rate=0.3,
                 use_time_aware=True, use_dnn=True):
        super(LSTM_SelfAttention_Model, self).__init__()
        
        # Control parameters
        self.use_time_aware = use_time_aware
        self.use_dnn = use_dnn
        
        # LSTM layers for workday and restday
        self.lstm_workday = nn.ModuleList(
            [nn.LSTM(input_size=1, hidden_size=hidden_size, batch_first=True) 
             for _ in range(input_size)]
        )
        self.lstm_restday = nn.ModuleList(
            [nn.LSTM(input_size=1, hidden_size=hidden_size, batch_first=True) 
             for _ in range(input_size)]
        )
        
        # Dropout layers
        self.lstm_dropout = nn.Dropout(lstm_dropout_rate)
        self.attention_dropout = nn.Dropout(attention_dropout_rate)
        
        # Self-attention layers
        self.self_attention_workday = nn.ModuleList(
            [nn.MultiheadAttention(embed_dim=hidden_size, num_heads=2) 
             for _ in range(input_size)]
        )
        self.self_attention_restday = nn.ModuleList(
            [nn.MultiheadAttention(embed_dim=hidden_size, num_heads=2) 
             for _ in range(input_size)]
        )
        
        # Time-aware module (optional)
        if self.use_time_aware:
            self.time_aware_module = nn.Sequential(
                nn.Linear(1, 16),
                nn.ReLU(),
                nn.Dropout(0.5),
                nn.Linear(16, 8),
                nn.ReLU(),
                nn.Dropout(0.5)
            )
        
        # DNN encoder for static features (optional)
        if self.use_dnn:
            self.dnn = nn.Sequential(
                nn.Linear(static_feature_size, 64),
                nn.ReLU(),
                nn.Dropout(0.5),
                nn.Linear(64, 32),
                nn.ReLU(),
                nn.Dropout(0.5)
            )
        
        # Calculate final fully connected layer input size
        fc_input_size = hidden_size * 2 * input_size  # Base LSTM features
        
        if self.use_dnn:
            fc_input_size += 32  # Add DNN feature dimension
        
        if self.use_time_aware:
            fc_input_size += 8  # Add time-aware feature dimension
            
        # Final classification layer
        self.fc = nn.Linear(fc_input_size, output_size)

    def self_attention(self, lstm_out, self_attention_layer):
        """Apply self-attention mechanism to LSTM output"""
        lstm_out = lstm_out.permute(1, 0, 2)
        attn_output, attn_weights = self_attention_layer(lstm_out, lstm_out, lstm_out)
        attn_output = attn_output.permute(1, 0, 2)
        weighted_sum = torch.mean(attn_output, dim=1)
        return weighted_sum, attn_weights
    
    def forward(self, workday_input, restday_input, static_input, time_gap):
        """Forward pass through the network"""
        # Store attention weights for analysis
        workday_attn_weights_list = []
        restday_attn_weights_list = []

        # Process workday branch
        workday_attended_list = []
        for i in range(workday_input.size(2)):
            workday_lstm_out, _ = self.lstm_workday[i](workday_input[:, :, i].unsqueeze(2))
            workday_lstm_out = self.lstm_dropout(workday_lstm_out)
            workday_attended, attn_weights = self.self_attention(
                workday_lstm_out, self.self_attention_workday[i]
            )
            workday_attended = self.attention_dropout(workday_attended)
            workday_attended_list.append(workday_attended)
            workday_attn_weights_list.append(attn_weights)

        # Process restday branch
        restday_attended_list = []
        for i in range(restday_input.size(2)):
            restday_lstm_out, _ = self.lstm_restday[i](restday_input[:, :, i].unsqueeze(2))
            restday_lstm_out = self.lstm_dropout(restday_lstm_out)
            restday_attended, attn_weights = self.self_attention(
                restday_lstm_out, self.self_attention_restday[i]
            )
            restday_attended = self.attention_dropout(restday_attended)
            restday_attended_list.append(restday_attended)
            restday_attn_weights_list.append(attn_weights)

        # Concatenate LSTM features
        combined_features = torch.cat(workday_attended_list + restday_attended_list, dim=1)
        combined_features = self.lstm_dropout(combined_features)

        # Process static features (if DNN is enabled)
        if self.use_dnn:
            static_features = self.dnn(static_input)
            combined_features = torch.cat((combined_features, static_features), dim=1)

        # Process time gap features (if time-aware is enabled)
        if self.use_time_aware:
            time_gap_features = self.time_aware_module(time_gap)
            combined_features = torch.cat((combined_features, time_gap_features), dim=1)

        # Final prediction
        output = self.fc(combined_features)

        return output, workday_attn_weights_list, restday_attn_weights_list