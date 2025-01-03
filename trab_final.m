%TRABALHO FINAL
%GRUPO 4 - IGOR MARTINS SILVA (00333069) E ARTUR RUIZ DE SOUZA (00334954)
close all;
clear;
clc;
%declaração das variáveis
mod_QPSK = 4; 
mod_QAM = 64;
multTamMensagem = 50 * 2;
N = 1296; %numero de bits da mensagem
n = N * multTamMensagem;  %tamanho da mensagem com código
R = 2/3; %razão de código
K = (N*R); %bits que são efeticamente a mensagem
k = R*n; %bits de mensagem de toda a transmissão
bits_pari = n - k; %bits de paridade
frame_size = 2300 * 8;
num_frames = k/K; % Número de quadros simulados por Eb/N0

Eb_N0_dB = -2:1:25; % Faixa de Eb/N0 em dB
Eb_N0_lin = 10 .^ (Eb_N0_dB / 10); % Faixa de Eb/N0 em linearizada

%Gerar informação
mensagem = randi(2,k,1)-1;
mensagemQAM = 0:(mod_QAM-1); %vetor para a consetelação do QAM
mensagemCod = zeros(n, 1);

%declaração das matrizes e vetores
ber_qpsk = zeros(3, length(Eb_N0_lin));
fer_qpsk = zeros(3, length(Eb_N0_lin));
ber_qam = zeros(3, length(Eb_N0_lin));
fer_qam = zeros(3, length(Eb_N0_lin));

%--------------------Matriz de paridade----------------------
 PariMatrix = qc_matrix_1296(N,R);
% PariMatrix = dvbs2ldpc(R);
%--------------------Codificação LDPC------------------------
%objeto codificador
ldpcEncoder = comm.LDPCEncoder(PariMatrix);
%objeto decodificador hard
ldpcDecoderHard = comm.LDPCDecoder('ParityCheckMatrix',PariMatrix, 'DecisionMethod','Hard decision'); 
%objeto decodificador soft
ldpcDecoderSoft = comm.LDPCDecoder('ParityCheckMatrix',PariMatrix, 'DecisionMethod','Soft decision', 'MaximumIterationCount', 50);

%codificação da mensagem em blocos de tamanho K=864 por conta das restrições da
%biblioteca COMM.ldpc
for j = 1: num_frames
    if j==1
        mensagemCod =  ldpcEncoder.step(mensagem((j-1)*K+1:j*K));
    else
        mensagemCod = cat(1, mensagemCod, ldpcEncoder.step(mensagem((j-1)*K+1:j*K)));
    end
end

%--------------------Modulação 64_qam------------------------
%Criação de uma constelação padrão 64-QAM para modulação da mensagem
qamMod = qammod(mensagemQAM, mod_QAM, 0 ,'Gray');
const = qamMod/ sqrt(mean(abs(qamMod).^2)); % Normalização para potência média unitária
QAMmod = comm.GeneralQAMModulator(const);

symbols_Cod = bi2de(reshape(mensagemCod, log2(mod_QAM), []).', 'left-msb');
qamModCod = qammod(symbols_Cod, mod_QAM, 0 ,'Gray');
qamModCod = qamModCod / sqrt(mean(abs(qamModCod).^2)); % Normalização para potência média unitária
% qamModCod = QAMmod.step(mensagemCod);
% qamModCod = qamModCod / sqrt(mean(abs(qamModCod).^2)); 

symbols = bi2de(reshape(mensagem, log2(mod_QAM), []).', 'left-msb');
qamMod = qammod(symbols, mod_QAM, 0 ,'Gray');
qamMod = qamMod / sqrt(mean(abs(qamMod).^2)); % Normalização para potência média unitária

% qamMod = QAMmod.step(mensagem);
% qamMod = qamMod / sqrt(mean(abs(qamMod).^2)); 

%----------------Objetos de Demodulação 64_qam-----------------
QAMdemod = comm.GeneralQAMDemodulator(const, 'BitOutput',true,'DecisionMethod','Hard decision');
QAMdemodHard = comm.GeneralQAMDemodulator(const, 'BitOutput',true,'DecisionMethod','Hard decision');
QAMdemodSoft = comm.GeneralQAMDemodulator(const, 'BitOutput',true,'DecisionMethod','Approximate log-likelihood ratio');

%--------------------Modulação QPSK----------------------------
qpskmod = comm.PSKModulator(mod_QPSK, 'BitInput',true);
qpskdemod = comm.PSKDemodulator(mod_QPSK, 'BitOutput',true,'DecisionMethod','Hard decision');
qpskdemodHard = comm.PSKDemodulator(mod_QPSK, 'BitOutput',true,'DecisionMethod','Hard decision');
qpskdemodSoft = comm.PSKDemodulator(mod_QPSK, 'BitOutput',true,'DecisionMethod','Approximate log-likelihood ratio');

qpsk_Mod = qpskmod.step(mensagem);
qpsk_Mod_Cod = qpskmod.step(mensagemCod);

%--------------------Cálculo BER para QPSK------------------
Eb = 1/log2(mod_QPSK); % Energia média para QPSK
NP = Eb ./ (Eb_N0_lin); %vetor de potências do ruído
NA = sqrt(NP); %vetor de amplitudes do ruído

EbCod = Eb/R; % Valores considerando a razão de código
NPCod = EbCod ./ (Eb_N0_lin);
NACod = sqrt(NPCod);
mensagemDemodDecodHard = zeros(k,1);
mensagemDemodDecodSoft = zeros(k,1);
for i = 1:length(Eb_N0_lin)
    NSemCod = NA(i)*complex(randn(length(qpsk_Mod), 1), randn(length(qpsk_Mod), 1))*sqrt(0.5); %vetor de rúido complexo com desvio padrão igual a uma posição do vetor NA
    NCod = NACod(i)*complex(randn(length(qpsk_Mod_Cod), 1), randn(length(qpsk_Mod_Cod), 1))*sqrt(0.5);
    rSemCod = qpsk_Mod + NSemCod; % vetor recebido sem código
    rCod = qpsk_Mod_Cod + NCod; % vetor recebido codificado

    mensagemDemod = qpskdemod.step(rSemCod);
    auxHard = 4-8.*qpskdemodHard.step(rCod);
    auxSoft = qpskdemodSoft.step(rCod);
    
    for j = 1: n/N
        if j==1
            mensagemDemodDecodHard =  ldpcDecoderHard.step(auxHard((j-1)*N+1:j*N));
            mensagemDemodDecodSoft = ldpcDecoderSoft.step(auxSoft((j-1)*N+1:j*N));
        else
            mensagemDemodDecodHard = cat(1, mensagemDemodDecodHard, ldpcDecoderHard.step(auxHard((j-1)*N+1:j*N)));
            mensagemDemodDecodSoft = cat(1, mensagemDemodDecodSoft, ldpcDecoderSoft.step(auxSoft((j-1)*N+1:j*N)));
        end
    end
    
    mensagemDemodDecodSoft = (sign(mensagemDemodDecodSoft)-1)/-2;
    

    ber_qpsk(1, i) = sum(mensagem ~= mensagemDemod) / k; % contagem de erros e cálculo do BER para QPSK sem codificação
    ber_qpsk(2, i) = sum(mensagem ~= mensagemDemodDecodHard) / k; % contagem de erros e cálculo do BER para QPSK codificação Hard
    ber_qpsk(3, i) = sum(mensagem ~= mensagemDemodDecodSoft) / k; % contagem de erros e cálculo do BER para QPSK codificação Soft
    % Calculo da probabilidade de um frame ter um ou mais erros
    fer_qpsk(1, i) = 1-((1-ber_qpsk(1,i))^frame_size);  
    fer_qpsk(2, i) = 1-((1-ber_qpsk(2,i))^frame_size);
    fer_qpsk(3, i) = 1-((1-ber_qpsk(3,i))^frame_size);
end

% %--------------------Cálculo BER para 64-QAM------------------        
EbQAM = mean(abs(const))/log2(mod_QAM); % Energia média para 64 - QAM
NPQAM = EbQAM ./ (Eb_N0_lin); %vetor de potências do ruido
NAQAM = sqrt(NPQAM); %vetor de amplitudes do ruído

EbCodQAM = EbQAM/R; % Valores considerando a razão de código
NPCodQAM = EbCodQAM ./ (Eb_N0_lin);
NACodQAM = sqrt(NPCodQAM);
mensagemDemodDecodHardQAM = zeros(k,1);
mensagemDemodDecodSoftQAM = zeros(k,1);
for i = 1:length(Eb_N0_lin)
    NSemCodQAM = NAQAM(i)*complex(randn(length(qamMod), 1), randn(length(qamMod), 1))*sqrt(0.5); %vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    NCodQAM = NACodQAM(i)*complex(randn(length(qamModCod), 1), randn(length(qamModCod), 1))*sqrt(0.5);
    rSemCodQAM = qamMod + NSemCodQAM; % vetor recebido
    rCodQAM = qamModCod + NCodQAM;
    
    mensagemDemodQAM = step(QAMdemod,rSemCodQAM);    
    auxHardQAM = 4-8.*QAMdemodHard.step(rCodQAM);
    auxSoftQAM = QAMdemodSoft.step(rCodQAM);
    
    for j = 1: n/N
        if j==1
            mensagemDemodDecodHardQAM =  ldpcDecoderHard.step(auxHardQAM((j-1)*N+1:j*N));
            mensagemDemodDecodSoftQAM = ldpcDecoderSoft.step(auxSoftQAM((j-1)*N+1:j*N));
        else
            mensagemDemodDecodHardQAM = cat(1, mensagemDemodDecodHardQAM, ldpcDecoderHard.step(auxHardQAM((j-1)*N+1:j*N)));
            mensagemDemodDecodSoftQAM = cat(1, mensagemDemodDecodSoftQAM, ldpcDecoderSoft.step(auxSoftQAM((j-1)*N+1:j*N)));
        end
    end
%   mensagemDemodDecodSoftQAM = (sign(mensagemDemodDecodSoftQAM)-1)/-2;
    mensagemDemodDecodSoftQAM = double(mensagemDemodDecodSoftQAM < 0);
    
    ber_qam(1, i) = sum(mensagem ~= mensagemDemodQAM) / k; % contagem de erros e calculo do BER para 64-QAM sem codificação
    ber_qam(2, i) = sum(mensagem ~= mensagemDemodDecodHardQAM) / k; % contagem de erros e calculo do BER para 64-QAM codificação Hard
    ber_qam(3, i) = sum(mensagem ~= mensagemDemodDecodSoftQAM) / k; % contagem de erros e calculo do BER para 64-QAM codificação Soft
    % Calculo da probabilidade de um frame ter um ou mais erros
    fer_qam(1, i) = 1-((1-ber_qam(1,i))^frame_size);
    fer_qam(2, i) = 1-((1-ber_qam(2,i))^frame_size);
    fer_qam(3, i) = 1-((1-ber_qam(3,i))^frame_size);
end


figure(1);
semilogy(Eb_N0_dB, ber_qpsk(1,:), 'r', 'LineWidth', 3); hold on;
semilogy(Eb_N0_dB, ber_qpsk(2,:), 'g', 'LineWidth', 3);
semilogy(Eb_N0_dB, ber_qpsk(3,:), 'b', 'LineWidth', 3);
semilogy(Eb_N0_dB, ber_qam(1,:), 'c', 'LineWidth', 3);
semilogy(Eb_N0_dB, ber_qam(2,:), 'k', 'LineWidth', 3); 
semilogy(Eb_N0_dB, ber_qam(3,:), 'm', 'LineWidth', 3); 
hold off;
grid on;
xlabel('Eb/N0 (dB)');
ylabel('BER');
legend('QPSK sem Cod', 'QPSK LDPC Hard', 'QPSK LDPC Soft', ...
       '64-QAM sem Cod', '64-QAM LDPC Hard', '64-QAM LDPC Soft');
title('Comparação BER de cada tipo de decodificação LDPC');

figure(2);
semilogy(Eb_N0_dB, fer_qpsk(1,:), 'r', 'LineWidth', 3); hold on;
semilogy(Eb_N0_dB, fer_qpsk(2,:), 'g', 'LineWidth', 4);
semilogy(Eb_N0_dB, fer_qpsk(3,:), 'b', 'LineWidth', 3);
hold off;
xlabel('Eb/N0 (dB)');
ylabel('FER');
legend('QPSK sem Cod', 'QPSK LDPC Hard', 'QPSK LDPC Soft');
title('Comparação Frame Error Rate QPSK');

figure(3);
semilogy(Eb_N0_dB, fer_qam(1,:), 'c', 'LineWidth', 3);hold on;
semilogy(Eb_N0_dB, fer_qam(2,:), 'k', 'LineWidth', 4); 
semilogy(Eb_N0_dB, fer_qam(3,:), 'm', 'LineWidth', 3); 
hold off;
xlabel('Eb/N0 (dB)');
ylabel('FER');
legend('64-QAM sem Cod', '64-QAM LDPC Hard', '64-QAM LDPC Soft');
title('Comparação Frame Error Rate QAM');


figure(4);
semilogy(Eb_N0_dB, fer_qpsk(2,:), 'g', 'LineWidth', 3);hold on;
semilogy(Eb_N0_dB, fer_qam(2,:), 'k', 'LineWidth', 4); 
hold off;
xlabel('Eb/N0 (dB)');
ylabel('FER');
legend('QPSK LDPC Hard', '64-QAM LDPC Hard');
title('Comparação Frame Error Rate HARD');
