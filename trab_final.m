%TRABALHO FINAL
%GRUPO 4 - IGOR MARTINS SILVA (00333069) E ARTUR RUIZ DE SOUZA (00334954)
close all;
clear;
clc;
%declara��o de vari�veis
mod_QPSK = 4; 
mod_QAM = 64;
n = 1296 * 50; %comprimento do c�digo
R = 2/3;
k = R*n; %bits de mensagem
bits_pari = n - k; %bits de paridade
frame_size = 2300 * 8;
Eb_N0_dB = -2:1:12; % Faixa de Eb/N0 em dB
Eb_N0_lin = 10 .^ (Eb_N0_dB / 10); % Faixa de Eb/N0 em linearizada
num_frames = 100; % N�mero de quadros simulados por Eb/N0

%Gerar informa��o
mensagem = randi(2,k,1)-1;
mensagemDemod = zeros(k, 1);
mensagemDemodDecod = zeros(k, 1);

%declara��o das matrizes e vetores
ber_qpsk = zeros(3, length(Eb_N0_lin));
fer_qpsk = zeros(3, length(Eb_N0_lin));
ber_qam = zeros(3, length(Eb_N0_lin));
fer_qam = zeros(3, length(Eb_N0_lin));

%--------------------Codifica��o LDPC------------------------
PariMatrix = dvbs2ldpc(R);
ldpcEncoder = comm.LDPCEncoder(PariMatrix);
ldpcDecoderHard = comm.LDPCDecoder('ParityCheckMatrix',PariMatrix, 'DecisionMethod','Hard decision');
ldpcDecoderSoft = comm.LDPCDecoder('ParityCheckMatrix',PariMatrix, 'DecisionMethod','Soft decision');
mensagemCod = ldpcEncoder.step(mensagem);

%--------------------Modula��o 64_qam------------------------
m1 = 0;  n1 = 0;
x = zeros(64,1);
for k1 = -7:2:7
    m1 = m1+1;
    for l = -7:2:7
        m1 = m1+1;
     x(m1,1) =  k1+1i*l;
   end 
   n1 = 0;
end

QAMmod = comm.GeneralQAMModulator(x);
QAMdemod = comm.GeneralQAMDemodulator(x);
QAMdemodHard = comm.GeneralQAMDemodulator(x, 'BitOutput',true,'DecisionMethod','Hard decision');
QAMdemodSoft = comm.GeneralQAMDemodulator(x, 'BitOutput',true,'DecisionMethod','Approximate log-likelihood ratio');

qamMod = QAMmod.step(mensagem);
qamModCod = QAMmod.step(mensagemCod);
% symbols = bi2de(reshape(mensagem, log2(mod_QAM), []).', 'left-msb');
% qamMod = qammod(symbols, mod_QAM, 0 ,'Gray');
% qamMod = qamMod / sqrt(mean(abs(qamMod).^2)); % Normaliza��o para pot�ncia m�dia unit�ria
% 
% 
% symbols_Cod = bi2de(reshape(mensagemCod, log2(mod_QAM), []).', 'left-msb');
% qamModCod = qammod(symbols_Cod, mod_QAM, 0 ,'Gray');
% qamModCod = qamModCod / sqrt(mean(abs(qamModCod).^2)); % Normaliza��o para pot�ncia m�dia unit�ria

% Plot da constela��o QPSK
%scatterplot(qamMod);
%title('Constela��o QPSK');

%se precisar fazer dar nos pontos 1,3,5 e 7 multiplicar por sqrt(42)

%--------------------Modula��o qpsk------------------------
% Agrupar bits e converter para inteiros (0 a 3)
% symbols_QPSK = bi2de(reshape(mensagem, log2(mod_QPSK), []).', 'left-msb');
% 
% symbols_QPSK_Cod = bi2de(reshape(mensagemCod, log2(mod_QPSK), []).', 'left-msb');

% Modula��o QPSK


% qpsk_Mod = pskmod(symbols_QPSK, mod_QPSK, pi/4); % pi/4 para constela��o usual
% 
% qpsk_Mod_Cod = pskmod(symbols_QPSK, mod_QPSK, pi/4); %modula��o dos bits codificados com LDPC


%scatterplot(qpsk_Mod);

qpskmod = comm.PSKModulator(mod_QPSK, 'BitInput',true);
qpskdemod = comm.PSKDemodulator(mod_QPSK, 'BitOutput',true,'DecisionMethod','Hard decision');
qpskdemodHard = comm.PSKDemodulator(mod_QPSK, 'BitOutput',true,'DecisionMethod','Hard decision');
qpskdemodSoft = comm.PSKDemodulator(mod_QPSK, 'BitOutput',true,'DecisionMethod','Approximate log-likelihood ratio');

qpsk_Mod = qpskmod.step(mensagem);
qpsk_Mod_Cod = qpskmod.step(mensagemCod);

%--------------------C�lculo BER para QPSK------------------
Eb = sqrt(2); % Energia m�dia para QPSK
NP = Eb ./ (Eb_N0_lin); %vetor de pot�ncias do ru�do
NA = sqrt(NP); %vetor de amplitudes do ru�do

EbCod = Eb/R; % Valores considerando a raz�o de c�digo
NPCod = EbCod ./ (Eb_N0_lin);
NACod = sqrt(NPCod);

for i = 1:length(Eb_N0_lin)
    NSemCod = NA(i)*complex(randn(length(qpsk_Mod), 1), randn(length(qpsk_Mod), 1))*sqrt(0.5); %vetor de ru�do complexo com desvio padr�o igual a uma posi��o do vetor NA
    NCod = NACod(i)*complex(randn(length(qpsk_Mod_Cod), 1), randn(length(qpsk_Mod_Cod), 1))*sqrt(0.5);
    rSemCod = qpsk_Mod + NSemCod; % vetor recebido
    rCod = qpsk_Mod_Cod + NCod;
    mensagemDemod = qpskdemod.step(rSemCod);
    mensagemDemodDecodHard = ldpcDecoderHard.step(qpskdemodHard.step(rCod));
    
    mensagemDemodDecodSoft = ldpcDecoderSoft.step(qpskdemodSoft.step(rCod))
    
    ber_qpsk(1, i) = sum(mensagem ~= mensagemDemod) / k; % contagem de erros e c�lculo do BER para QPSK sem codifica��o
    ber_qpsk(2, i) = sum(mensagem ~= mensagemDemodDecodHard) / k; % contagem de erros e c�lculo do BER QPSK com codifica��o Hard
    ber_qpsk(3, i) = sum(mensagem ~= mensagemDemodDecodSoft) / k; % contagem de erros e c�lculo do BER QPSK com codifica��o Soft
end

figure(1);
semilogy(Eb_N0_dB, ber_qpsk(1,:), 'r', Eb_N0_dB, ber_qpsk(2,:), 'g', Eb_N0_dB, ber_qpsk(3,:), 'b');
xlabel('Eb/N0 (dB)');
ylabel('BER');
legend('QPSK sem Cod', 'QPSK LDPC Hard', 'QPSK LDPC Soft');
% 
% 
% 
% 
% %--------------------C�lculo BER para 64-QAM------------------
% Eb = sqrt(2); % Energia m�dia para 64 - QAM (???????)  
% NP = Eb ./ (Eb_N0_lin); %vetor de pot�ncias do ru�do
% NA = sqrt(NP); %vetor de amplitudes do ru�do
% 
% EbCod = Eb/R; % Valores considerando a raz�o de c�digo
% NPCod = EbCod ./ (Eb_N0_lin);
% NACod = sqrt(NPCod);
% 
% for i = 1:length(Eb_N0_lin)
%     NSemCod = NA(i)*complex(randn(length(qamMod), 1), randn(length(qamMod), 1))*sqrt(0.5); %vetor de ru�do complexo com desvio padr�o igual a uma posi��o do vetor NA
%     NCod = NACod(i)*complex(randn(length(qamModCod), 1), randn(length(qamModCod), 1))*sqrt(0.5)*3/4;
%     rSemCod = qamMod + NSemCod; % vetor recebido
%     rCod = qamModCod + NCod;
%     mensagemDemod = QAMdemod.step(rSemCod);
%     mensagemDemodDecodHard = ldpcDecoderHard.step(QAMdemodHard.step(rCod));
%     mensagemDemodDecodSoft = ldpcDecoderSoft.step(QAMdemodSoft.step(rCod));
%     
%     ber_qam(1, i) = sum(mensagem ~= mensagemDemod) / k; % contagem de erros e c�lculo do BER para 64-QAM sem codifica��o
%     ber_qam(2, i) = sum(mensagem ~= mensagemDemodDecodHard) / k; % contagem de erros e c�lculo do BER 64-QAM com codifica��o Hard
%     ber_qam(3, i) = sum(mensagem ~= mensagemDemodDecodSoft) / k; % contagem de erros e c�lculo do BER 64-QAM com codifica��o Soft
% end
% 
% figure(2);
% semilogy(Eb_N0_dB, ber_qam(1,:), 'r', Eb_N0_dB, ber_qam(2,:), 'g', Eb_N0_dB, ber_qam(3,:), 'b');
% xlabel('Eb/N0 (dB)');
% ylabel('BER');
% legend('64-QAM sem Cod', '64-QAM LDPC Hard', '64-QAM LDPC Soft');

% %--------------------------------------LDPC encoder------------------------
% 
% 
% % Gerar bits de mensagem aleat�rios
% message_bits = randi([0 1], k, 1);
% 
% % Codificar
% encoded_bits = ldpc_encode(message_bits, H, n, k, bits_pari);
% 
% % Verificar dimens�es
% %disp(['Bits codificados: ', num2str(length(encoded_bits))]); % Deve ser igual a n
% 
% 
% % Mensagem aleat�ria
% %input_bits = randi([0 1], 120, 1); % 120 bits (20 s�mbolos para 64-QAM)
% 
% % Modula��o
% qam_symbols = manual_qammod(message_bits,mod_QAM);
% 
% % Visualizar constela��o
% % Exibir apenas os pontos da constela��o
% figure(2);
% scatter(real(qam_symbols), imag(qam_symbols), 'filled');
% xlabel('Eixo I (In-Phase)');
% ylabel('Eixo Q (Quadrature)');
% title('Constela��o 64-QAM');
% grid on;
% axis([-8 8 -8 8]); % Ajustar os limites para visualiza��o correta
% 












