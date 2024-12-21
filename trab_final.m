%TRABALHO FINAL
%GRUPO 4 - IGOR MARTINS SILVA (00333069) E ARTUR RUIZ DE SOUZA (00334954)
close all;
clear;
clc;
%declaração de variáveis
mod_QPSK = 4; 
mod_QAM = 64;
n = 1296 * 50; %comprimento do código
R = 2/3;
k = R*n; %bits de mensagem
bits_pari = n - k; %bits de paridade
frame_size = 2300 * 8;
Eb_N0_dB = -2:1:12; % Faixa de Eb/N0 em dB
Eb_N0_lin = 10 .^ (Eb_N0_dB / 10); % Faixa de Eb/N0 em linearizada
num_frames = 100; % Número de quadros simulados por Eb/N0

%Gerar informação
mensagem = randi(2,k,1)-1;
mensagemQAM = 0:(mod_QAM-1);
mensagemDemod = zeros(k, 1);
mensagemDemodDecod = zeros(k, 1);

%declaração das matrizes e vetores
ber_qpsk = zeros(3, length(Eb_N0_lin));
fer_qpsk = zeros(3, length(Eb_N0_lin));
ber_qam = zeros(3, length(Eb_N0_lin));
fer_qam = zeros(3, length(Eb_N0_lin));

%--------------------Codificação LDPC------------------------
PariMatrix = dvbs2ldpc(R);
ldpcEncoder = comm.LDPCEncoder(PariMatrix);
ldpcDecoderHard = comm.LDPCDecoder('ParityCheckMatrix',PariMatrix, 'DecisionMethod','Hard decision');
ldpcDecoderSoft = comm.LDPCDecoder('ParityCheckMatrix',PariMatrix, 'DecisionMethod','Hard decision');
mensagemCod = ldpcEncoder.step(mensagem);

%--------------------Modulação 64_qam------------------------
% symbols = bi2de(reshape(mensagem, log2(mod_QAM), []).', 'left-msb');
qamMod = qammod(mensagemQAM, mod_QAM, 0 ,'Gray');
const = qamMod/ sqrt(mean(abs(qamMod).^2)); % Normalização para potência média unitária


symbols_Cod = bi2de(reshape(mensagemCod, log2(mod_QAM), []).', 'left-msb');
qamModCod = qammod(symbols_Cod, mod_QAM, 0 ,'Gray');
qamModCod = qamModCod / sqrt(mean(abs(qamModCod).^2)); % Normalização para potência média unitária

symbols = bi2de(reshape(mensagem, log2(mod_QAM), []).', 'left-msb');
qamMod = qammod(symbols, mod_QAM, 0 ,'Gray');
qamMod = qamMod / sqrt(mean(abs(qamMod).^2)); % Normalização para potência média unitária

QAMmod = comm.GeneralQAMModulator(const);
QAMdemod = comm.GeneralQAMDemodulator(const, 'BitOutput',true,'DecisionMethod','Hard decision');
QAMdemodHard = comm.GeneralQAMDemodulator(const, 'BitOutput',true,'DecisionMethod','Hard decision');
QAMdemodSoft = comm.GeneralQAMDemodulator(const, 'BitOutput',true,'DecisionMethod','Approximate log-likelihood ratio');

% qamMod = QAMmod.step(mensagem);
% qamModCod = QAMmod.step(mensagemCod);
% Plot da constelação QPSK
% scatterplot(qamMod);
% title('Constelação QPSK');

%se precisar fazer dar nos pontos 1,3,5 e 7 multiplicar por sqrt(42)

%--------------------Modulação qpsk------------------------
% Agrupar bits e converter para inteiros (0 a 3)
% symbols_QPSK = bi2de(reshape(mensagem, log2(mod_QPSK), []).', 'left-msb');
% 
% symbols_QPSK_Cod = bi2de(reshape(mensagemCod, log2(mod_QPSK), []).', 'left-msb');

% Modulação QPSK


% qpsk_Mod = pskmod(symbols_QPSK, mod_QPSK, pi/4); % pi/4 para constelação usual
% 
% qpsk_Mod_Cod = pskmod(symbols_QPSK, mod_QPSK, pi/4); %modulação dos bits codificados com LDPC


%scatterplot(qpsk_Mod);

qpskmod = comm.PSKModulator(mod_QPSK, 'BitInput',true);
qpskdemod = comm.PSKDemodulator(mod_QPSK, 'BitOutput',true,'DecisionMethod','Hard decision');
qpskdemodHard = comm.PSKDemodulator(mod_QPSK, 'BitOutput',true,'DecisionMethod','Hard decision');
qpskdemodSoft = comm.PSKDemodulator(mod_QPSK, 'BitOutput',true,'DecisionMethod','Approximate log-likelihood ratio');

qpsk_Mod = qpskmod.step(mensagem);
qpsk_Mod_Cod = qpskmod.step(mensagemCod);

%--------------------Cálculo BER para QPSK------------------
Eb = 1; % Energia média para QPSK
NP = Eb ./ (Eb_N0_lin); %vetor de potências do ruído
NA = sqrt(NP); %vetor de amplitudes do ruído

EbCod = Eb*R; % Valores considerando a razão de código
NPCod = EbCod ./ (Eb_N0_lin);
NACod = sqrt(NPCod);

for i = 1:length(Eb_N0_lin)
    NSemCod = NA(i)*complex(randn(length(qpsk_Mod), 1), randn(length(qpsk_Mod), 1))*sqrt(0.5); %vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    NCod = NACod(i)*complex(randn(length(qpsk_Mod_Cod), 1), randn(length(qpsk_Mod_Cod), 1))*sqrt(0.5);
    rSemCod = qpsk_Mod + NSemCod; % vetor recebido
    rCod = qpsk_Mod_Cod + NCod;
    mensagemDemod = qpskdemod.step(rSemCod);
    auxHard = 4-8.*qpskdemodHard.step(rCod);
    mensagemDemodDecodHard = ldpcDecoderHard.step(auxHard);
    
    mensagemDemodDecodSoft = ldpcDecoderSoft.step(qpskdemodSoft.step(rCod));
    
    ber_qpsk(1, i) = sum(mensagem ~= mensagemDemod) / k; % contagem de erros e cálculo do BER para QPSK sem codificação
    ber_qpsk(2, i) = sum(mensagem ~= mensagemDemodDecodHard) / k; % contagem de erros e cálculo do BER QPSK com codificação Hard
    ber_qpsk(3, i) = sum(mensagem ~= mensagemDemodDecodSoft) / k; % contagem de erros e cálculo do BER QPSK com codificação Soft
end




% %--------------------Cálculo BER para 64-QAM------------------
Eb = mean(abs(const)); % Energia média para 64 - QAM
NP = Eb ./ (Eb_N0_lin); %vetor de potências do ruído
NA = sqrt(NP); %vetor de amplitudes do ruído

EbCod = Eb*R; % Valores considerando a razão de código
NPCod = EbCod ./ (Eb_N0_lin);
NACod = sqrt(NPCod);

for i = 1:length(Eb_N0_lin)
    NSemCod = NA(i)*complex(randn(length(qamMod), 1), randn(length(qamMod), 1))*sqrt(0.5); %vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    NCod = NACod(i)*complex(randn(length(qamModCod), 1), randn(length(qamModCod), 1))*sqrt(0.5)*3/4;
    rSemCod = qamMod + NSemCod; % vetor recebido
    rCod = qamModCod + NCod;
    
    mensagemDemod = QAMdemod.step(rSemCod);
    
    auxHard = 4-8.*QAMdemodHard.step(rCod);
    auxSoft = QAMdemodSoft.step(rCod);
    
    mensagemDemodDecodHard = ldpcDecoderHard.step(auxHard);
    mensagemDemodDecodSoft = ldpcDecoderSoft.step(auxSoft);
  
    
    ber_qam(1, i) = sum(mensagem ~= mensagemDemod) / k; % contagem de erros e cálculo do BER para 64-QAM sem codificação
    ber_qam(2, i) = sum(mensagem ~= mensagemDemodDecodHard) / k; % contagem de erros e cálculo do BER 64-QAM com codificação Hard
    ber_qam(3, i) = sum(mensagem ~= mensagemDemodDecodSoft) / k; % contagem de erros e cálculo do BER 64-QAM com codificação Soft
end

% figure(1);
% semilogy(Eb_N0_dB, ber_qam(1,:), 'r', Eb_N0_dB, ber_qam(2,:), 'g', Eb_N0_dB, ber_qam(3,:), 'b');
% xlabel('Eb/N0 (dB)');
% ylabel('BER');
% legend('64-QAM sem Cod', '64-QAM LDPC Hard', '64-QAM LDPC Soft');





figure(1);
semilogy(Eb_N0_dB, ber_qpsk(1,:), 'r', Eb_N0_dB, ber_qpsk(2,:), 'g', Eb_N0_dB, ber_qpsk(3,:), 'b', Eb_N0_dB, ber_qam(1,:), 'y', Eb_N0_dB, ber_qam(2,:), 'black', Eb_N0_dB, ber_qam(3,:), 'magenta');
xlabel('Eb/N0 (dB)');
ylabel('BER');
legend('QPSK sem Cod', 'QPSK LDPC Hard', 'QPSK LDPC Soft','64-QAM sem Cod', '64-QAM LDPC Hard', '64-QAM LDPC Soft');










