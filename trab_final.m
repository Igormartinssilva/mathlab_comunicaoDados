%TRABALHO FINAL
%GRUPO 4 - IGOR MARTINS SILVA (00333069) E ARTUR RUIZ DE SOUZA (00334954)
close all;
clear;
clc;
%declaração de variáveis
mod_QPSK = 4; 
mod_QAM = 64;
n = 64800; %comprimento do código
num_b = n; 
R = 2/3;
k = R*n; %bits de mensagem
bits_pari = n - k; %bits de paridade
frame_size = 2300 * 8;
Eb_N0_dB = -2:1:12; % Faixa de Eb/N0 em dB
Eb_N0_lin = 10 .^ (Eb_N0_dB / 10); % Faixa de Eb/N0 em linearizada
num_frames = 100; % Número de quadros simulados por Eb/N0

mensagem = randi(2,k,1)-1;
mensagemDemod = zeros(k, 1);
mensagemDemodDecod = zeros(k, 1);

%declaração das matrizes e vetores
ber_qpsk = zeros(3, length(Eb_N0_lin));
fer_qpsk = zeros(3, length(Eb_N0_lin));
ber_qam = zeros(3, length(Eb_N0_lin));
fer_qam = zeros(3, length(Eb_N0_lin));


% QPSK ---------------------------------------------------------

% Criar matriz H na forma [P | I]
P = randi([0 1], bits_pari, k); % Parte aleatória (paridade)
I = eye(bits_pari);             % Matriz identidade
H = [P I];              % Combinação de P e I


%Variáveis
Eb = sqrt(2); % Energia média para quadratura
NP = Eb ./ (Eb_N0_lin); %vetor de potências do ruído
NA = sqrt(NP); %vetor de amplitudes do ruído
ldpcEncoder = comm.LDPCEncoder(dvbs2ldpc(R));
ldpcDecoder = comm.LDPCDecoder('ParityCheckMatrix',dvbs2ldpc(R), 'DecisionMethod','Hard decision');
qpskmod = comm.QPSKModulator;
qpskdemodHard = comm.QPSKDemodulator;
qpskdemodSoft = comm.QPSKDemodulator('DecisionMethod', 'Log-likelihood ratio');

mensagemCod = ldpcEncoder.step(mensagem);
SModCod = qpskmod.step(mensagemCod); 
SMod = qpskmod.step(mensagem); %Sinal resultante, modulado e complexo

for i = 1:length(Eb_N0_lin)
    NSemCod = NA(i)*complex(randn(length(SMod), 1), randn(length(SMod), 1))*sqrt(0.5); %vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    NCod = NA(i)*complex(randn(length(SModCod), 1), randn(length(SModCod), 1))*sqrt(0.5);
    rSemCod = SMod + NSemCod; % vetor recebido
    rCod = SModCod + NCod;
    
    mensagemDemod = qpskdemodHard.step(rSemCod);
    mensagemDemodDecodHard = ldpcDecoder.step(qpskdemodHard.step(rCod));
    mensagemDemodDecodSoft = ldpcDecoder.step(qpskdemodSoft.step(rCod));
    
    ber_qpsk(1, i) = sum(mensagem ~= mensagemDemod) / k; % contagem de erros e cálculo do BER para QPSK sem codificação
    ber_qpsk(2, i) = sum(mensagem ~= mensagemDemodDecodHard) / k; % contagem de erros e cálculo do BER QPSK com codificação Hard
    ber_qpsk(3, i) = sum(mensagem ~= mensagemDemodDecodSoft) / k; % contagem de erros e cálculo do BER QPSK com codificação Soft
end

figure(1);
semilogy(Eb_N0_dB, ber_qpsk(1,:), 'r', Eb_N0_dB, ber_qpsk(2,:), 'g', Eb_N0_dB, ber_qpsk(3,:), 'b');
xlabel('Eb/N0 (dB)');
ylabel('BER');
legend('QPSK sem Cod', 'QPSK LDPC Hard', 'QPSK LDPC Soft');

%--------------------------------------LDPC encoder------------------------


% Gerar bits de mensagem aleatórios
message_bits = randi([0 1], k, 1);

% Codificar
encoded_bits = ldpc_encode(message_bits, H, n, k, bits_pari);

% Verificar dimensões
%disp(['Bits codificados: ', num2str(length(encoded_bits))]); % Deve ser igual a n


% Mensagem aleatória
%input_bits = randi([0 1], 120, 1); % 120 bits (20 símbolos para 64-QAM)

% Modulação
qam_symbols = manual_qammod(message_bits,mod_QAM);

% Visualizar constelação
% Exibir apenas os pontos da constelação
figure(2);
scatter(real(qam_symbols), imag(qam_symbols), 'filled');
xlabel('Eixo I (In-Phase)');
ylabel('Eixo Q (Quadrature)');
title('Constelação 64-QAM');
grid on;
axis([-8 8 -8 8]); % Ajustar os limites para visualização correta













