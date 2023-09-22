%% Relay commuication amplify and forward method

clear all;
close all;
clc;

Ps = 1;     % Source power
Pr = Ps/2;  % Relay power
N = 10^6;   % Number of bits

eta = 4;    % Pathloss exponent

d_sr = 100;  % Distance from source to relay
d_rd = 100;  % Distance from relay to destination
d_sd = 200;  % Distance from source to destination

snr_db = 0:2:40;                % SNR range in dB
snr = (10^-3)*db2pow(snr_db);   % SNR range in linear scale

data = randi([0 1],1,N);
bpsk_data = 2*binary_data - 1;

% Channel gains
h_sr = sqrt(1000*(d_sr^-eta))*(1/sqrt(2))*(randn(1,N)+1i*randn(1,N));
h_rd = sqrt(1000*(d_rd^-eta))*(1/sqrt(2))*(randn(1,N)+1i*randn(1,N));
h_sd = sqrt(1000*(d_sd^-eta))*(1/sqrt(2))*(randn(1,N)+1i*randn(1,N));

% Noises
n_sr = (randn(1,N) + 1i*randn(1,N))/sqrt(2);
n_rd = (randn(1,N) + 1i*randn(1,N))/sqrt(2);
n_sd = (randn(1,N) + 1i*randn(1,N))/sqrt(2);

