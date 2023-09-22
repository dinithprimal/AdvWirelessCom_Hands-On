clc;
close all
clear all;

% Parameters
numBits = 100000; % Number of bits
SNR_dB = 0:2:20;  % SNR range in dB

% Convert SNR from dB to linear scale
SNR_linear = 10.^(SNR_dB / 10);

% Relay amplification gain
G = 2;

% Initialize BER vector
BER = zeros(size(SNR_dB));

% Generate random bits
bits = randi([0 1], 1, numBits);

% BPSK modulation
symbols = 2 * bits - 1;

for snrIdx = 1:length(SNR_dB)
    
    
    
    % Transmit through source-relay-destination
    h_sr = (randn(1, numBits) + 1i * randn(1, numBits)) / sqrt(2);
    h_rd = (randn(1, numBits) + 1i * randn(1, numBits)) / sqrt(2);
    
    % Add AWGN at source
    Es = 1; % Symbol energy for BPSK
    noiseVar_source = Es / (2 * SNR_linear(snrIdx));
    noise_source = sqrt(noiseVar_source) * randn(size(symbols));
    received_source = symbols + noise_source;
    
    % Amplify the signal at the relay
    received_relay = G * received_source .* h_sr;
    
    % Add AWGN at relay
    noiseVar_relay = Es / (2 * SNR_linear(snrIdx));
    noise_relay = sqrt(noiseVar_relay) * randn(size(received_relay));
    received_relay = received_relay + noise_relay;
    
    % Amplify the signal at the destination
    received_destination = received_relay .* h_rd;
    
    % Add AWGN at destination
    noiseVar_destination = Es / (2 * SNR_linear(snrIdx));
    noise_destination = sqrt(noiseVar_destination) * randn(size(received_destination));
    received_destination = received_destination + noise_destination;
    
    % LLR calculation
    LLR = 2 * received_destination / noiseVar_destination;
    
    % LLR-based soft decision demodulation
    decodedBits = LLR < 0;
    
    % Calculate BER
    errors = sum(decodedBits ~= bits);
    BER(snrIdx) = errors / numBits;
end

% Plot BER vs SNR
semilogy(SNR_dB, BER, 'o-');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('AF Relay Channel with BPSK Modulation and LLR-based Soft Decision Demodulation');
grid on;
