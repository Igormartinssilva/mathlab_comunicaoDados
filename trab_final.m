%TRABALHO FINAL
%GRUPO 4 - IGOR MARTINS SILVA (00333069) E ARTUR RUIZ DE SOUZA (00334954)
close all;
clear;
clc;
%declara��o de vari�veis
mod_QPSK = 4; 
mod_QAM = 64;
n = 64800; %comprimento do c�digo
num_b = n; 
R = 2/3;
k = R*n; %bits de mensagem
bits_pari = n - k; %bits de paridade
frame_size = 2300 * 8;
Eb_N0_dB = -2:1:12; % Faixa de Eb/N0 em dB
Eb_N0_lin = 10 .^ (Eb_N0_dB / 10); % Faixa de Eb/N0 em linearizada
num_frames = 100; % N�mero de quadros simulados por Eb/N0

mensagem = randi(2,k,1)-1;
mensagemDemod = zeros(k, 1);
mensagemDemodDecod = zeros(k, 1);

%declara��o das matrizes e vetores
ber_qpsk = zeros(3, length(Eb_N0_lin));
fer_qpsk = zeros(3, length(Eb_N0_lin));
ber_qam = zeros(3, length(Eb_N0_lin));
fer_qam = zeros(3, length(Eb_N0_lin));


% QPSK ---------------------------------------------------------

% Criar matriz H na forma [P | I]
P = randi([0 1], bits_pari, k); % Parte aleat�ria (paridade)
I = eye(bits_pari);             % Matriz identidade
H = [P I];              % Combina��o de P e I


%Vari�veis
Eb = sqrt(2); % Energia m�dia para quadratura
NP = Eb ./ (Eb_N0_lin); %vetor de pot�ncias do ru�do
NA = sqrt(NP); %vetor de amplitudes do ru�do
ldpcEncoder = comm.LDPCEncoder(dvbs2ldpc(R));
ldpcDecoder = comm.LDPCDecoder('ParityCheckMatrix',dvbs2ldpc(R), 'DecisionMethod','Hard decision');
qpskmod = comm.QPSKModulator;
qpskdemodHard = comm.QPSKDemodulator;
qpskdemodSoft = comm.QPSKDemodulator('DecisionMethod', 'Log-likelihood ratio');

mensagemCod = ldpcEncoder.step(mensagem);
SModCod = qpskmod.step(mensagemCod); 
SMod = qpskmod.step(mensagem); %Sinal resultante, modulado e complexo

for i = 1:length(Eb_N0_lin)
    NSemCod = NA(i)*complex(randn(length(SMod), 1), randn(length(SMod), 1))*sqrt(0.5); %vetor de ru�do complexo com desvio padr�o igual a uma posi��o do vetor NA
    NCod = NA(i)*complex(randn(length(SModCod), 1), randn(length(SModCod), 1))*sqrt(0.5);
    rSemCod = SMod + NSemCod; % vetor recebido
    rCod = SModCod + NCod;
    
    mensagemDemod = qpskdemodHard.step(rSemCod);
    mensagemDemodDecodHard = ldpcDecoder.step(qpskdemodHard.step(rCod));
    mensagemDemodDecodSoft = ldpcDecoder.step(qpskdemodSoft.step(rCod));
    
    ber_qpsk(1, i) = sum(mensagem ~= mensagemDemod) / k; % contagem de erros e c�lculo do BER para QPSK sem codifica��o
    ber_qpsk(2, i) = sum(mensagem ~= mensagemDemodDecodHard) / k; % contagem de erros e c�lculo do BER QPSK com codifica��o Hard
    ber_qpsk(3, i) = sum(mensagem ~= mensagemDemodDecodSoft) / k; % contagem de erros e c�lculo do BER QPSK com codifica��o Soft
end

figure(1);
semilogy(Eb_N0_dB, ber_qpsk(1,:), 'r', Eb_N0_dB, ber_qpsk(2,:), 'g', Eb_N0_dB, ber_qpsk(3,:), 'b');
xlabel('Eb/N0 (dB)');
ylabel('BER');
legend('QPSK sem Cod', 'QPSK LDPC Hard', 'QPSK LDPC Soft');

%--------------------------------------LDPC encoder------------------------


% Gerar bits de mensagem aleat�rios
message_bits = randi([0 1], k, 1);

% Codificar
encoded_bits = ldpc_encode(message_bits, H, n, k, bits_pari);

% Verificar dimens�es
%disp(['Bits codificados: ', num2str(length(encoded_bits))]); % Deve ser igual a n


% Mensagem aleat�ria
%input_bits = randi([0 1], 120, 1); % 120 bits (20 s�mbolos para 64-QAM)

% Modula��o
qam_symbols = manual_qammod(message_bits,mod_QAM);

% Visualizar constela��o
% Exibir apenas os pontos da constela��o
figure(2);
scatter(real(qam_symbols), imag(qam_symbols), 'filled');
xlabel('Eixo I (In-Phase)');
ylabel('Eixo Q (Quadrature)');
title('Constela��o 64-QAM');
grid on;
axis([-8 8 -8 8]); % Ajustar os limites para visualiza��o correta













