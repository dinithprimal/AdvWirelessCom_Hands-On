clc;
clear all;
close all;

%% System model parameters

N = 10^6;

d = 500;        % Distance from source to destination
d_S_R = d/2;    % Distance from source to Relay
d_R_D = d/2;    % Distance from Relay to destination

eta = 4;    % Path Loss Exponent
k = 1;

h_D = sqrt(d^-eta)*k;%(randn(1,1) + 1i*randn(1,1))/sqrt(2);        %source to destination channel
h_S_R = sqrt(d_S_R^-eta)*k;%(randn(1,1) + 1i*randn(1,1))/sqrt(2);  %source to relay channel
h_S_R = sqrt(d_S_R^-eta)*k;%(randn(1,1) + 1i*randn(1,1))/sqrt(2);  %relay to destination channel

%Channel gains
g_D = (abs(h_D)).^2;
g_S_R = (abs(h_S_R)).^2;
g_R_D = (abs(h_S_R)).^2;

%% Simulation

SNR = 0:40;         %SNR in dB
snr = db2pow(SNR);  %SNR linear

% SNR = P/(sigma^2), Transmit SNR
% where P = Total Transmit Power
% sigma = noise power

%snr_S_D = zeros(1,length(snr));
%snr_S_R_D = zeros(1,length(snr));

for u = 1:length(snr)
    snr_S_D(u) = pow2db(snr(u)*g_D);
    snr_S_R_D(u) = pow2db(0.5*snr(u)*g_S_R + 0.5*snr(u)*g_R_D);
    %gain_R_D(u) = 10*log10(0.5*snr(u)*g_R_D);
end

plot(SNR,snr_S_D,'linewidth',2); hold on; grid on;
plot(SNR,snr_S_R_D,'linewidth',2);

xlabel('Transmitted SNR in dB'); ylabel('Received SNR in dB');
legend('Direct Communication','Relay Communication');


