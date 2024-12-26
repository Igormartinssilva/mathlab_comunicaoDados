%TRABALHO FINAL
%GRUPO 4 - IGOR MARTINS SILVA (00333069) E ARTUR RUIZ DE SOUZA (00334954)
close all;
clear;
clc;
%declaração de variáveis
mod_QPSK = 4; 
mod_QAM = 64;
tamanhoMensagem = 3;
n = 1296 * 50 * tamanhoMensagem; %comprimento do código
R = 2/3;
k = R*n; %bits de mensagem
bits_pari = n - k; %bits de paridade
frame_size = 2300 * 8;
Eb_N0_dB = -2:1:12; % Faixa de Eb/N0 em dB
Eb_N0_lin = 10 .^ (Eb_N0_dB / 10); % Faixa de Eb/N0 em linearizada
num_frames = n/frame_size; % Número de quadros simulados por Eb/N0

%Gerar informação
mensagem = randi(2,k,1)-1;
mensagemQAM = 0:(mod_QAM-1);
mensagemDemod = zeros(k, 1);
mensagemCod = zeros(n, 1);
mensagemDemodDecod = zeros(k, 1);

%declaração das matrizes e vetores
ber_qpsk = zeros(3, length(Eb_N0_lin));
fer_qpsk = zeros(3, length(Eb_N0_lin));
ber_qam = zeros(3, length(Eb_N0_lin));
fer_qam = zeros(3, length(Eb_N0_lin));

%--------------------Matriz de paridade----------------------
PariMatrix = ldpcMatrizDeParidade(1296*50,R);

%----Fazer Gráfico de densidade, bom gráfico, talvez usar no futuro-----

% figure(1);
% imagesc(PariMatrix(1:1000, 1:1000)); % Exibe as primeiras 50 linhas e colunas
% colormap(gray);         % Colormap em tons de cinza
% colorbar;               % Adiciona uma barra de cores
% title('Trecho da matriz H (50x50)');
% xlabel('Colunas');
% ylabel('Linhas');
%--------------------Codificação LDPC------------------------

%teste = dvbs2ldpc(R);
ldpcEncoder = comm.LDPCEncoder(PariMatrix);
ldpcDecoderHard = comm.LDPCDecoder('ParityCheckMatrix',PariMatrix, 'DecisionMethod','Hard decision');
ldpcDecoderSoft = comm.LDPCDecoder('ParityCheckMatrix',PariMatrix, 'DecisionMethod','Soft decision');

for j = 1: k/43200
    if j==1
         mensagemCod =  ldpcEncoder.step(mensagem((j-1)*43200+1:j*43200));
    else
        mensagemCod = cat(1, mensagemCod, ldpcEncoder.step(mensagem((j-1)*43200+1:j*43200)));
    end
end

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
mensagemDemodDecodHard = zeros(k,1);
mensagemDemodDecodSoft = zeros(k,1);
for i = 1:length(Eb_N0_lin)
    NSemCod = NA(i)*complex(randn(length(qpsk_Mod), 1), randn(length(qpsk_Mod), 1))*sqrt(0.5); %vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    NCod = NACod(i)*complex(randn(length(qpsk_Mod_Cod), 1), randn(length(qpsk_Mod_Cod), 1))*sqrt(0.5);
    rSemCod = qpsk_Mod + NSemCod; % vetor recebido
    rCod = qpsk_Mod_Cod + NCod;
    mensagemDemod = qpskdemod.step(rSemCod);
    
    auxHard = 4-8.*qpskdemodHard.step(rCod);
    auxSoft = qpskdemodSoft.step(rCod);
    
    for j = 1: n/64800
        if j==1
            mensagemDemodDecodHard =  ldpcDecoderHard.step(auxHard((j-1)*64800+1:j*64800));
            mensagemDemodDecodSoft = ldpcDecoderSoft.step(auxSoft((j-1)*64800+1:j*64800));
        else
            mensagemDemodDecodHard = cat(1, mensagemDemodDecodHard, ldpcDecoderHard.step(auxHard((j-1)*64800+1:j*64800)));
            mensagemDemodDecodSoft = cat(1, mensagemDemodDecodSoft, ldpcDecoderSoft.step(auxSoft((j-1)*64800+1:j*64800)));
        end
    end
    
    mensagemDemodDecodSoft = (sign(mensagemDemodDecodSoft)-1)/-2;
    index = floor(num_frames);
    for j = 1: index
        if sum(mensagem((j-1)*index+1:j*index)~= mensagemDemod((j-1)*index+1:j*index) > 0)
            fer_qpsk(1, i) = fer_qpsk(1, i) + 1;
        end
        
        if sum(mensagem((j-1)*index+1:j*index)~= mensagemDemodDecodHard((j-1)*index+1:j*index) > 0)
            fer_qpsk(2, i) = fer_qpsk(1, i) + 1;
        end
        
        if sum(mensagem((j-1)*index+1:j*index)~= mensagemDemodDecodSoft((j-1)*index+1:j*index) > 0)
            fer_qpsk(3, i) = fer_qpsk(1, i) + 1;
        end
    end
    
    fer_qpsk(1, i) = fer_qpsk(1, i)/num_frames;
    fer_qpsk(2, i) = fer_qpsk(1, i)/num_frames;
    fer_qpsk(3, i) = fer_qpsk(1, i)/num_frames;
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
    
    for j = 1: n/64800
        if j==1
            mensagemDemodDecodHard =  ldpcDecoderHard.step(auxHard((j-1)*64800+1:j*64800));
            mensagemDemodDecodSoft = ldpcDecoderSoft.step(auxSoft((j-1)*64800+1:j*64800));
        else
            mensagemDemodDecodHard = cat(1, mensagemDemodDecodHard, ldpcDecoderHard.step(auxHard((j-1)*64800+1:j*64800)));
            mensagemDemodDecodSoft = cat(1, mensagemDemodDecodSoft, ldpcDecoderSoft.step(auxSoft((j-1)*64800+1:j*64800)));
        end
    end
  
    mensagemDemodDecodSoft = (sign(mensagemDemodDecodSoft)-1)/-2;
    index = floor(num_frames);
    for j = 1: index
        if sum(mensagem((j-1)*index+1:j*index)~= mensagemDemod((j-1)*index+1:j*index) > 0)
            fer_qam(1, i) = fer_qam(1, i) + 1;
        end
        
        if sum(mensagem((j-1)*index+1:j*index)~= mensagemDemodDecodHard((j-1)*index+1:j*index) > 0)
            fer_qam(2, i) = fer_qam(1, i) + 1;
        end
        
        if sum(mensagem((j-1)*index+1:j*index)~= mensagemDemodDecodSoft((j-1)*index+1:j*index) > 0)
            fer_qam(3, i) = fer_qam(1, i) + 1;
        end
    end
    
    fer_qam(1, i) = fer_qam(1, i)/num_frames;
    fer_qam(2, i) = fer_qam(1, i)/num_frames;
    fer_qam(3, i) = fer_qam(1, i)/num_frames;
    
    ber_qam(1, i) = sum(mensagem ~= mensagemDemod) / k; % contagem de erros e cálculo do BER para 64-QAM sem codificação
    ber_qam(2, i) = sum(mensagem ~= mensagemDemodDecodHard) / k; % contagem de erros e cálculo do BER 64-QAM com codificação Hard
    ber_qam(3, i) = sum(mensagem ~= mensagemDemodDecodSoft) / k; % contagem de erros e cálculo do BER 64-QAM com codificação Soft
end


figure(1);
semilogy(Eb_N0_dB, ber_qpsk(1,:), 'r', 'LineWidth', 2); hold on;
semilogy(Eb_N0_dB, ber_qpsk(2,:), 'g', 'LineWidth', 2);
semilogy(Eb_N0_dB, ber_qpsk(3,:), 'b', 'LineWidth', 2);
semilogy(Eb_N0_dB, ber_qam(1,:), 'c', 'LineWidth', 2);
semilogy(Eb_N0_dB, ber_qam(2,:), 'k', 'LineWidth', 2); 
semilogy(Eb_N0_dB, ber_qam(3,:), 'm', 'LineWidth', 2); 
hold off;
xlabel('Eb/N0 (dB)');
ylabel('BER');
legend('QPSK sem Cod', 'QPSK LDPC Hard', 'QPSK LDPC Soft', ...
       '64-QAM sem Cod', '64-QAM LDPC Hard', '64-QAM LDPC Soft');

figure(2);
semilogy(Eb_N0_dB, fer_qpsk(1,:), 'r', 'LineWidth', 2); hold on;
semilogy(Eb_N0_dB, fer_qpsk(2,:), 'g', 'LineWidth', 2);
semilogy(Eb_N0_dB, fer_qpsk(3,:), 'b', 'LineWidth', 2);
semilogy(Eb_N0_dB, fer_qam(1,:), 'c', 'LineWidth', 2);
semilogy(Eb_N0_dB, fer_qam(2,:), 'k', 'LineWidth', 2); 
semilogy(Eb_N0_dB, fer_qam(3,:), 'm', 'LineWidth', 2); 
hold off;
xlabel('Eb/N0 (dB)');
ylabel('FER');
legend('QPSK sem Cod', 'QPSK LDPC Hard', 'QPSK LDPC Soft', ...
       '64-QAM sem Cod', '64-QAM LDPC Hard', '64-QAM LDPC Soft');
