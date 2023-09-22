clc;
clear all
close all;

N = 10^5; % Number of bits

d_SR = 0.5; 
d_RD = 0.5;
d_SD = d_SR + d_RD;

Pt = 40; % Transmit power in dBm
pt = (10^-3)*db2pow(Pt); % Transmit power in linear scale

Nr = 1; % Flat fading

% SNR slace
SNR = 0:2:40;			%SNR (dBm)
snr = (10^-3)*db2pow(SNR);	%SNR (linear scale)


eta = 4;

% Rayliegh fading channel
h_SR = sqrt(d_SR^-eta)*RayleighFading(Nr);
h_RD = sqrt(d_RD^-eta)*RayleighFading(Nr);
h_SD = sqrt(d_SD^-eta)*RayleighFading(Nr);


BW = 10^6;			%Bandwidth = 1 MHz
No = -174 + 10*log10(BW);	%Noise power (dBm)			%
no = (10^-3)*db2pow(No);	%Noise power (linear scale) 

%Generate noise samples for the three paths
n_SR = sqrt(no)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
n_RD = sqrt(no)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
n_SD = sqrt(no)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);

x0 = randi([0 1],N,1);
% BPSK modulation
x = 2 * x0 - 1;

y_SR = sqrt(pt/2)*x.*h_SR + n_SR;
eq_y_SR = y_SR./h_SR;


for u = 1:length(snr)
    
    noisePow = pt./snr(u);
    
    beta = sqrt(pt./(((abs(h_SR).^2))+noisePow));
    
    
    
    y_RD = sqrt(pt/2)*(beta.*eq_y_SR).*h_RD + n_RD;
    eq_y_RD = y_RD./h_RD;
    
    y_SD = sqrt(pt)*x.*h_SD + n_SD;
    eq_y_SD = y_SD./h_SD;
    
    LLR_SD = (2*sqrt(pt(u))*y_SD.*h_SD)/var(n_SD);
    %LLR_SD = (2*y_SD)./n_SD;
    eq_LLR_SD = LLR_SD./h_SD;
    
    LLR_T = (2*y_RD.*beta.*(h_SR.*h_RD))./(var(n_RD)*(1+(pt(u).*(abs(h_RD).^2).*(beta.^2)))) + LLR_SD;
    eq_LLR_T = LLR_T./h_RD;
    
    x_RD = zeros(N,1);% dummy array to store the far user decoded bits
    x_RD(eq_y_RD>0) = 1;
    
    x_SD = zeros(N,1);% dummy array to store the far user decoded bits
    x_SD(eq_y_SD>0) = 1;
    
    x_LLR_SD = zeros(N,1);% dummy array
    x_LLR_SD(eq_LLR_SD>0) = 1;
    
    x_LLR_T = zeros(N,1);% dummy array to store the far user decoded bits
    x_LLR_T(eq_LLR_T>0) = 1;
    
    ber_RD(u) = biterr(x0,x_RD)/N;
    ber_SD(u) = biterr(x0,x_SD)/N;
    ber_LLR_SD(u) = biterr(x0,x_LLR_SD)/N;
    ber_LLR_T(u) = biterr(x0,x_LLR_T)/N;
    
end

semilogy(Pt,ber_RD,'b'),hold on,
semilogy(Pt,ber_SD,'r')
semilogy(Pt,ber_LLR_SD,'c')
semilogy(Pt,ber_LLR_T,'m')





