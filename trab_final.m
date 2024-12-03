%TRABALHO FINAL
%GRUPO 4 - IGOR MARTINS SILVA (00333069) E ARTUR RUIZ DE SOUZA (00334954)
close all;
clear;

%delcaração de variáveis
mod_QPSK = 4; 
mod_QUAM = 64;
n = 1296;
R = 2/3;
frame_size = 2300 * 8;
Eb_N0_dB = -2:1:12; % Faixa de Eb/N0 em dB
Eb_N0_lin = 10 .^ (Eb_N0_dB / 10); % Faixa de Eb/N0 em linearizada
num_frames = 100; % Número de quadros simulados por Eb/N0


NP = Eb ./ (Eb_N0_lin); %vetor de potências do ruído
NA = sqrt(NP); %vetor de amplitudes do ruído


%declaração das matrizes e vetores
ber_qpsk = zeros(3, length(Eb_N0_dB));
fer_qpsk = zeros(3, length(Eb_N0_dB));
ber_qam = zeros(3, length(Eb_N0_dB));
fer_qam = zeros(3, length(Eb_N0_dB));




























