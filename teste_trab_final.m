% TRABALHO FINAL
% GRUPO 4 - IGOR MARTINS SILVA (00333069) E ARTUR RUIZ DE SOUZA (00334954)
close all;
clear;
clc;

% Declaracao de variaveis
mod_QPSK = 4;
mod_QAM = 64;
num_bits_QPSK = log2(mod_QPSK); % numero de bits por simbolo
num_bits_QAM = log2(mod_QAM); % numero de bits por simbolo
n = 1296; % numero de bits de cada mensagem
r = 2/3; % taxa de correcao
k = n * r; % bits que sao realmente a mensagem
bit_pari = n - k; % bits que fazem a correcao da mensagem
frame_size = 2300 * 8; % numero de bits em um frame (2300 bytes)
num_frames = 100; % numero de frames/blocos
num_bits = frame_size * num_frames; % numero de bits
num_simbols_QPSK = num_bits / num_bits_QPSK; % numero de simbolos
num_simbols_QAM = num_bits / num_bits_QAM; % numero de simbolos

% Calculo da energia de bit
Eb_N0_db = -2:1:12; % intervalo do eb
Eb_N0_lin = 10 .^ (Eb_N0_db / 10); % Eb linearizado

% Geracao das mensagens
mensagem = randi([0 1], num_bits, 1); % gera a mensagem 
mensagem_Cod = []; % mensagem codificada

% Geracao dos vetores BER
BER_NoCoding = zeros(1, length(Eb_N0_lin));
BER_Soft = zeros(1, length(Eb_N0_lin));
BER_Hard = zeros(1, length(Eb_N0_lin));

% Codificaï¿½ï¿½o LDPC
% p = [
%         39, 31, 22, 43, -1, 40,  4, -1, 11, -1, -1, 50, -1, -1, -1,  6,  1,  0, -1, -1, -1, -1, -1, -1;
%         25, 52, 41,  2,  6, -1, 14, -1, 34, -1, -1, -1, 24, -1, 37, -1, -1,  0,  0, -1, -1, -1, -1, -1;
%         43, 31, 29,  0, 21, -1, 28, -1, -1,  2, -1, -1,  7, -1, 17, -1, -1, -1,  0,  0, -1, -1, -1, -1;
%         20, 33, 48, -1,  4, 13, -1, 26, -1, -1, 22, -1, -1, 46, 42, -1, -1, -1, -1,  0,  0, -1, -1, -1;
%         45,  7, 18, 51, 12, 25, -1, -1, -1, 50, -1, -1,  5, -1, -1, -1,  0, -1, -1, -1,  0,  0, -1, -1;
%         35, 40, 32, 16,  5, -1, -1, 18, -1, -1, 43, 51, -1, 32, -1, -1, -1, -1, -1, -1, -1,  0,  0, -1;
%          9, 24, 13, 22, 28, -1, -1, 37, -1, -1, 25, -1, -1, 52, -1, 13, -1, -1, -1, -1, -1, -1,  0,  0;
%         32, 22,  4, 21, 16, -1, -1, -1, 27, 28, -1, 38, -1, -1, -1,  8,  1, -1, -1, -1, -1, -1, -1,  0;
%     ];
% blockSize=54;
% H = ldpcQuasiCyclicMatrix(blockSize,p);
% % spy(H);
H = logical(qc_matrix_1296(n, r)); % Matriz de paridade esparsa e quasi-ciclica
%H = dvbs2ldpc(r);
%  H = ldpcMatrizDeParidade(n,r);
spy(H);
% Criaï¿½ï¿½o dos objetos ldpc
ldpcEncoder = comm.LDPCEncoder(H); % Codificador
ldpcDecoderHard = comm.LDPCDecoder(H); % Decodificador hard
ldpcDecoderSoft = comm.LDPCDecoder(H, 'DecisionMethod', 'Soft decision'); % Decodificador soft

% Codificar a mensagem em blocos
for i = 1:num_frames % Para cada frame codifica a mensagem
    
    bloco = mensagem((i-1)*k+1:i*k);    % Seleciona o bloco atual da mensagem de tamanho k=894
    bloco_codificado = step(ldpcEncoder, bloco);    % Codifica o bloco usando o LDPCEncoder   
    mensagem_Cod = [mensagem_Cod; bloco_codificado];    % Concatena o bloco codificado ï¿½ mensagem codificada
end

qpskmod = comm.QPSKModulator(mod_QPSK,'BitInput', true);
qpskDemod= comm.QPSKDemodulator(mod_QPSK,'BitOutput', true, 'DecisionMethod', 'hard decision');
qpskDemodHard= comm.QPSKDemodulator(mod_QPSK,'BitOutput', true, 'DecisionMethod', 'hard decision');
qpskDemodSoft= comm.QPSKDemodulator(mod_QPSK,'BitOutput', true, 'DecisionMethod', 'Log-likelihood ratio');

mensagem_mod_qpsk= qpskmod.step(mensagem);
mensagem_mod_qpsk_cod = qpskmod.step(mensagem_Cod);

%--------------calculo do BER-------------
Eb = 1 / log2(mod_QPSK); % Energia de bit
NP = Eb ./ Eb_N0_lin; % Potencia do ruido
NA = sqrt(NP / 2); % Desvio padrao do ruido

EbCod = Eb / r; % Energia de bit
NPCod = EbCod ./ Eb_N0_lin; % Potencia do ruido
NACod = sqrt(NPCod / 2); % Desvio padrao do ruido



for idx = 1:length(Eb_N0_db)
    % Calcular a SNR para o canal AWGN
    SNR = Eb_N0_db(idx) + 10*log10(num_bits_QPSK);

    % Adiciona ruï¿½do para simulaï¿½ï¿½o do canal (nï¿½o codificado)
    mensagem_noisy = awgn(double(mensagem_mod_qpsk), SNR, 'measured');
    mensagem_dec_no_coding = mensagem_noisy;
  
    mensageDemod = qpskDemod.step(mensagem_dec_no_coding);
    BER_NoCoding(idx) = biterr(mensagem, mensageDemod) / num_bits;

    % Adiciona ruï¿½do para simulaï¿½ï¿½o do canal (codificado)
    mensagem_Cod_noisy = awgn(double(mensagem_mod_qpsk_cod), SNR, 'measured');

    % Decodificar a mensagem codificada
    mensagem_Dec_Hard = zeros(num_bits, 1); % mensagem decodificada hard
    mensagem_Demod_Soft = zeros(num_bits, 1); % mensagem decodificada soft

    auxHard = 4-8 .* qpskDemodHard.step(mensagem_Cod_noisy);
    auxSoft = qpskDemodSoft.step(mensagem_Cod_noisy);


    for i = 1:num_frames % Para cada frame decodifica a mensagem
        % Seleciona o bloco codificado
        bloco_codificado_noisy = auxHard((i-1)*n+1:i*n);

        % Decodifica o bloco usando o LDPCDecoder (hard decision)
        bloco_dec_hard = step(ldpcDecoderHard, bloco_codificado_noisy);
        mensagem_Dec_Hard((i-1)*k+1:i*k) = bloco_dec_hard;
        
        bloco_codificado_noisy_Soft = auxSoft((i-1)*n+1:i*n);
        % Ajusta para soft decision e decodifica
         
         bloco_noisy_soft =  bloco_codificado_noisy_Soft ;
        bloco_dec_soft = step(ldpcDecoderSoft, bloco_codificado_noisy_Soft);
        mensagem_Dec_Soft((i-1)*k+1:i*k) = bloco_dec_soft;
        mensagem_Dec_Soft = (sign(mensagem_Dec_Soft)-1) /-1;
    end

mensagem_Dec_Soft = mensagem_Dec_Soft(:); %faz a transposta da matriz
    % Verificaï¿½ï¿½o dos resultados
BER_Hard(idx) = biterr(mensagem(1:k*num_frames), mensagem_Dec_Hard(1:k*num_frames)) / (k*num_frames);
BER_Soft(idx) = biterr(mensagem(1:k*num_frames), mensagem_Dec_Soft(1:k*num_frames)) / (k*num_frames);
end

% Analisar a saída do decodificador soft
% figure;
% hist(bloco_dec_soft);
% title('Histograma da saída do decodificador soft');
% figure;
% hist(mensagem_Dec_Hard);
% title('Histograma da saída do decodificador hard');

BER_Hard
BER_Soft

% size(mensagem);
% size(mensagem_Dec_Hard);
% Verificaï¿½ï¿½o dos resultados
% [numErrorsHard, ratioHard] = biterr(mensagem, mensagem_Dec_Hard);
% [numErrorsSoft, ratioSoft] = biterr(mensagem, mensagem_Dec_Soft);
% [numErrorsUnCode, ratioUnCode] = biterr(mensagem/ num_bits, BER_NoCoding);
% for i = 1:num_frames
%     if i==1
%             erroHard(i) =  biterr(mensagem((i-1)*k+1:i*k), mensagem_Dec_Hard((i-1)*k+1:i*k));
%             erroSoft(i) =  biterr(mensagem((i-1)*k+1:i*k), mensagem_Dec_Soft((i-1)*k+1:i*k));
% 
%     else
%             erroHard(i) = erroHard(i-1) + biterr(mensagem((i-1)*k+1:i*k), mensagem_Dec_Hard((i-1)*k+1:i*k));
%             erroSoft(i) = erroSoft(i-1) + biterr(mensagem((i-1)*k+1:i*k), mensagem_Dec_Soft((i-1)*k+1:i*k));
% 
%     end
% end
% disp(['Numero de erros de bit (Hard Decision): ', num2str(erroHard)]);
% disp(['Numero de erros de bit (Soft Decision): ', num2str(erroSoft)]);
% 
% disp(['Numero de erros de bit (Hard Decision): ', num2str(numErrorsHard)]);
% disp(['Taxa de erros de bit (Hard Decision): ', num2str(ratioHard)]);
% disp(['Numero de erros de bit (Soft Decision): ', num2str(numErrorsSoft)]);
% disp(['Taxa de erros de bit (Soft Decision): ', num2str(ratioSoft)]);
% % disp(['Numero de erros de bit (Soft Decision): ', num2str(numErrorsUnCode)]);
% % disp(['Taxa de erros de bit (Soft Decision): ', num2str(ratioUnCode)]);



% Plotar os resultados
figure(2);
semilogy(Eb_N0_db, BER_NoCoding, '-r', 'DisplayName', 'Nï¿½o Codificado','LineWidth', 3);
hold on;
semilogy(Eb_N0_db, BER_Hard, '-b', 'DisplayName', 'Hard Decision','LineWidth', 3);
semilogy(Eb_N0_db, BER_Soft, '-g', 'DisplayName', 'Soft Decision','LineWidth', 3);
hold off;
grid on;
xlabel('Eb/N0 (dB)');
ylabel('BER');
legend('show');
title('Comparaï¿½ï¿½o de BER para Codificaï¿½ï¿½o LDPC');