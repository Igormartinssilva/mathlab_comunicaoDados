%TRABALHO FINAL
%GRUPO 4 - IGOR MARTINS SILVA (00333069) E ARTUR RUIZ DE SOUZA (00334954)
close all;
clear;
clc;
%declara��o das vari�veis
mod_QPSK = 4; 
mod_QAM = 64;
multTamMensagem = 50 * 2;
N = 1296; %numero de bits da mensagem
n = N * multTamMensagem;  %tamanho da mensagem com c�digo
R = 2/3; %raz�o de c�digo
K = (N*R); %bits que s�o efeticamente a mensagem
k = R*n; %bits de mensagem de toda a transmiss�o
bits_pari = n - k; %bits de paridade
frame_size = 2300 * 8;
num_frames = k/K; % N�mero de quadros simulados por Eb/N0

Eb_N0_dB = -2:1:25; % Faixa de Eb/N0 em dB
Eb_N0_lin = 10 .^ (Eb_N0_dB / 10); % Faixa de Eb/N0 em linearizada

%Gerar informa��o
mensagem = randi(2,k,1)-1;
mensagemQAM = 0:(mod_QAM-1); %vetor para a consetela��o do QAM
mensagemCod = zeros(n, 1);

%declara��o das matrizes e vetores
ber_qpsk = zeros(3, length(Eb_N0_lin));
fer_qpsk = zeros(3, length(Eb_N0_lin));
ber_qam = zeros(3, length(Eb_N0_lin));
fer_qam = zeros(3, length(Eb_N0_lin));

%--------------------Matriz de paridade----------------------
 PariMatrix = qc_matrix_1296(N,R);
% PariMatrix = dvbs2ldpc(R);
%--------------------Codifica��o LDPC------------------------
%objeto codificador
ldpcEncoder = comm.LDPCEncoder(PariMatrix);
%objeto decodificador hard
ldpcDecoderHard = comm.LDPCDecoder('ParityCheckMatrix',PariMatrix, 'DecisionMethod','Hard decision'); 
%objeto decodificador soft
ldpcDecoderSoft = comm.LDPCDecoder('ParityCheckMatrix',PariMatrix, 'DecisionMethod','Soft decision', 'MaximumIterationCount', 50);

%codifica��o da mensagem em blocos de tamanho K=864 por conta das restri��es da
%biblioteca COMM.ldpc
for j = 1: num_frames
    if j==1
        mensagemCod =  ldpcEncoder.step(mensagem((j-1)*K+1:j*K));
    else
        mensagemCod = cat(1, mensagemCod, ldpcEncoder.step(mensagem((j-1)*K+1:j*K)));
    end
end

%--------------------Modula��o 64_qam------------------------
%Cria��o de uma constela��o padr�o 64-QAM para modula��o da mensagem
qamMod = qammod(mensagemQAM, mod_QAM, 0 ,'Gray');
const = qamMod/ sqrt(mean(abs(qamMod).^2)); % Normaliza��o para pot�ncia m�dia unit�ria
QAMmod = comm.GeneralQAMModulator(const);

symbols_Cod = bi2de(reshape(mensagemCod, log2(mod_QAM), []).', 'left-msb');
qamModCod = qammod(symbols_Cod, mod_QAM, 0 ,'Gray');
qamModCod = qamModCod / sqrt(mean(abs(qamModCod).^2)); % Normaliza��o para pot�ncia m�dia unit�ria
% qamModCod = QAMmod.step(mensagemCod);
% qamModCod = qamModCod / sqrt(mean(abs(qamModCod).^2)); 

symbols = bi2de(reshape(mensagem, log2(mod_QAM), []).', 'left-msb');
qamMod = qammod(symbols, mod_QAM, 0 ,'Gray');
qamMod = qamMod / sqrt(mean(abs(qamMod).^2)); % Normaliza��o para pot�ncia m�dia unit�ria

% qamMod = QAMmod.step(mensagem);
% qamMod = qamMod / sqrt(mean(abs(qamMod).^2)); 

%----------------Objetos de Demodula��o 64_qam-----------------
QAMdemod = comm.GeneralQAMDemodulator(const, 'BitOutput',true,'DecisionMethod','Hard decision');
QAMdemodHard = comm.GeneralQAMDemodulator(const, 'BitOutput',true,'DecisionMethod','Hard decision');
QAMdemodSoft = comm.GeneralQAMDemodulator(const, 'BitOutput',true,'DecisionMethod','Approximate log-likelihood ratio');

%--------------------Modula��o QPSK----------------------------
qpskmod = comm.PSKModulator(mod_QPSK, 'BitInput',true);
qpskdemod = comm.PSKDemodulator(mod_QPSK, 'BitOutput',true,'DecisionMethod','Hard decision');
qpskdemodHard = comm.PSKDemodulator(mod_QPSK, 'BitOutput',true,'DecisionMethod','Hard decision');
qpskdemodSoft = comm.PSKDemodulator(mod_QPSK, 'BitOutput',true,'DecisionMethod','Approximate log-likelihood ratio');

qpsk_Mod = qpskmod.step(mensagem);
qpsk_Mod_Cod = qpskmod.step(mensagemCod);

%--------------------C�lculo BER para QPSK------------------
Eb = 1/log2(mod_QPSK); % Energia m�dia para QPSK
NP = Eb ./ (Eb_N0_lin); %vetor de pot�ncias do ru�do
NA = sqrt(NP); %vetor de amplitudes do ru�do

EbCod = Eb/R; % Valores considerando a raz�o de c�digo
NPCod = EbCod ./ (Eb_N0_lin);
NACod = sqrt(NPCod);
mensagemDemodDecodHard = zeros(k,1);
mensagemDemodDecodSoft = zeros(k,1);
for i = 1:length(Eb_N0_lin)
    NSemCod = NA(i)*complex(randn(length(qpsk_Mod), 1), randn(length(qpsk_Mod), 1))*sqrt(0.5); %vetor de r�ido complexo com desvio padr�o igual a uma posi��o do vetor NA
    NCod = NACod(i)*complex(randn(length(qpsk_Mod_Cod), 1), randn(length(qpsk_Mod_Cod), 1))*sqrt(0.5);
    rSemCod = qpsk_Mod + NSemCod; % vetor recebido sem c�digo
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
    

    ber_qpsk(1, i) = sum(mensagem ~= mensagemDemod) / k; % contagem de erros e c�lculo do BER para QPSK sem codifica��o
    ber_qpsk(2, i) = sum(mensagem ~= mensagemDemodDecodHard) / k; % contagem de erros e c�lculo do BER para QPSK codifica��o Hard
    ber_qpsk(3, i) = sum(mensagem ~= mensagemDemodDecodSoft) / k; % contagem de erros e c�lculo do BER para QPSK codifica��o Soft
    % Calculo da probabilidade de um frame ter um ou mais erros
    fer_qpsk(1, i) = 1-((1-ber_qpsk(1,i))^frame_size);  
    fer_qpsk(2, i) = 1-((1-ber_qpsk(2,i))^frame_size);
    fer_qpsk(3, i) = 1-((1-ber_qpsk(3,i))^frame_size);
end

% %--------------------C�lculo BER para 64-QAM------------------        
EbQAM = mean(abs(const))/log2(mod_QAM); % Energia m�dia para 64 - QAM
NPQAM = EbQAM ./ (Eb_N0_lin); %vetor de pot�ncias do ruido
NAQAM = sqrt(NPQAM); %vetor de amplitudes do ru�do

EbCodQAM = EbQAM/R; % Valores considerando a raz�o de c�digo
NPCodQAM = EbCodQAM ./ (Eb_N0_lin);
NACodQAM = sqrt(NPCodQAM);
mensagemDemodDecodHardQAM = zeros(k,1);
mensagemDemodDecodSoftQAM = zeros(k,1);
for i = 1:length(Eb_N0_lin)
    NSemCodQAM = NAQAM(i)*complex(randn(length(qamMod), 1), randn(length(qamMod), 1))*sqrt(0.5); %vetor de ru�do complexo com desvio padr�o igual a uma posi��o do vetor NA
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
    
    ber_qam(1, i) = sum(mensagem ~= mensagemDemodQAM) / k; % contagem de erros e calculo do BER para 64-QAM sem codifica��o
    ber_qam(2, i) = sum(mensagem ~= mensagemDemodDecodHardQAM) / k; % contagem de erros e calculo do BER para 64-QAM codifica��o Hard
    ber_qam(3, i) = sum(mensagem ~= mensagemDemodDecodSoftQAM) / k; % contagem de erros e calculo do BER para 64-QAM codifica��o Soft
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
title('Compara��o BER de cada tipo de decodifica��o LDPC');

figure(2);
semilogy(Eb_N0_dB, fer_qpsk(1,:), 'r', 'LineWidth', 3); hold on;
semilogy(Eb_N0_dB, fer_qpsk(2,:), 'g', 'LineWidth', 4);
semilogy(Eb_N0_dB, fer_qpsk(3,:), 'b', 'LineWidth', 3);
hold off;
xlabel('Eb/N0 (dB)');
ylabel('FER');
legend('QPSK sem Cod', 'QPSK LDPC Hard', 'QPSK LDPC Soft');
title('Compara��o Frame Error Rate QPSK');

figure(3);
semilogy(Eb_N0_dB, fer_qam(1,:), 'c', 'LineWidth', 3);hold on;
semilogy(Eb_N0_dB, fer_qam(2,:), 'k', 'LineWidth', 4); 
semilogy(Eb_N0_dB, fer_qam(3,:), 'm', 'LineWidth', 3); 
hold off;
xlabel('Eb/N0 (dB)');
ylabel('FER');
legend('64-QAM sem Cod', '64-QAM LDPC Hard', '64-QAM LDPC Soft');
title('Compara��o Frame Error Rate QAM');


figure(4);
semilogy(Eb_N0_dB, fer_qpsk(2,:), 'g', 'LineWidth', 3);hold on;
semilogy(Eb_N0_dB, fer_qam(2,:), 'k', 'LineWidth', 4); 
hold off;
xlabel('Eb/N0 (dB)');
ylabel('FER');
legend('QPSK LDPC Hard', '64-QAM LDPC Hard');
title('Compara��o Frame Error Rate HARD');
