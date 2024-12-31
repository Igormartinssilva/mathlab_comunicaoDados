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
H = qc_matrix_1296(n, r); % Matriz de paridade esparsa e quasi-ciclica
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










num_blocos = floor(num_bits / k); % Total de blocos necessï¿½rios
%decodificação sem ruido
for i = 1:num_frames
    % Seleciona o bloco codificado
    bloco_codificado = mensagem_Cod((i-1)*n+1 : i*n);

    % Decodifica usando Hard Decision
    bloco_decodificado_hard = step(ldpcDecoderHard,bloco_codificado);

    % Decodifica usando Soft Decision
    bloco_decodificado_soft = step(ldpcDecoderSoft,bloco_codificado);

    % Concatena os blocos decodificados
    mensagem_Dec_Hard((i-1)*k+1:i*k) = bloco_decodificado_hard;
    mensagem_Dec_Soft((i-1)*k+1:i*k) = bloco_decodificado_soft;
end
size(mensagem_Dec_Hard);
size(mensagem_Dec_Soft);
size(mensagem);
% Verifica erros
erroHard = sum(mensagem(1:length(mensagem_Dec_Hard(:))) ~= mensagem_Dec_Hard(:));
erroSoft = sum(mensagem(1:length(mensagem_Dec_Soft(:))) ~= mensagem_Dec_Soft(:));

disp(['Erros de bits no Hard Decision sem ruido: ', num2str(erroHard)]);
disp(['Erros de bits no Soft Decision sem ruido: ', num2str(erroSoft)]);






for idx = 1:length(Eb_N0_db)
    % Calcular a SNR para o canal AWGN
    SNR = Eb_N0_db(idx) + 10*log10(num_bits_QPSK);

    % Adiciona ruï¿½do para simulaï¿½ï¿½o do canal (nï¿½o codificado)
    mensagem_noisy = awgn(double(mensagem), SNR, 'measured');
    mensagem_dec_no_coding = mensagem_noisy > 0.5;
    BER_NoCoding(idx) = biterr(mensagem, mensagem_dec_no_coding) / num_bits;

    % Adiciona ruï¿½do para simulaï¿½ï¿½o do canal (codificado)
    mensagem_Cod_noisy = awgn(double(mensagem_Cod), SNR, 'measured');

    % Decodificar a mensagem codificada
    mensagem_Dec_Hard = zeros(num_bits, 1); % mensagem decodificada hard
    mensagem_Dec_Soft = zeros(num_bits, 1); % mensagem decodificada soft

    for i = 1:num_frames % Para cada frame decodifica a mensagem
        % Seleciona o bloco codificado
        bloco_codificado_noisy = mensagem_Cod_noisy((i-1)*n+1:i*n);
        
        % Decodifica o bloco usando o LDPCDecoder (hard decision)
        bloco_dec_hard = step(ldpcDecoderHard, bloco_codificado_noisy);
        mensagem_Dec_Hard((i-1)*k+1:i*k) = bloco_dec_hard;
        
        % Ajusta para soft decision e decodifica
        bloco_noisy_soft = 2 * bloco_codificado_noisy - 1; % Ajusta para valores -1 e 1
        bloco_dec_soft = step(ldpcDecoderSoft, bloco_noisy_soft);
        mensagem_Dec_Soft((i-1)*k+1:i*k) = bloco_dec_soft;
    end

    % Conversï¿½o para valores binï¿½rios (0 e 1) se necessï¿½rio
    %mensagem_Dec_Hard = mensagem_Dec_Hard > 0.5;
    mensagem_Dec_Soft = mensagem_Dec_Soft > 0.5;


    % Verificaï¿½ï¿½o dos resultados
    BER_Hard(idx) = biterr(mensagem, mensagem_Dec_Hard) / num_bits;
    BER_Soft(idx) = biterr(mensagem, mensagem_Dec_Soft) / num_bits;
end

%mensagem_Dec_Soft = mensagem_Dec_Soft > 0.5;

size(mensagem);
size(mensagem_Dec_Hard);
% Verificaï¿½ï¿½o dos resultados
[numErrorsHard, ratioHard] = biterr(mensagem, mensagem_Dec_Hard);
[numErrorsSoft, ratioSoft] = biterr(mensagem, mensagem_Dec_Soft);
% [numErrorsUnCode, ratioUnCode] = biterr(mensagem/ num_bits, BER_NoCoding);
for i = 1:num_frames
    if i==1
            erroHard(i) =  biterr(mensagem((i-1)*k+1:i*k), mensagem_Dec_Hard((i-1)*k+1:i*k));
            erroSoft(i) =  biterr(mensagem((i-1)*k+1:i*k), mensagem_Dec_Soft((i-1)*k+1:i*k));

    else
            erroHard(i) = erroHard(i-1) + biterr(mensagem((i-1)*k+1:i*k), mensagem_Dec_Hard((i-1)*k+1:i*k));
            erroSoft(i) = erroSoft(i-1) + biterr(mensagem((i-1)*k+1:i*k), mensagem_Dec_Soft((i-1)*k+1:i*k));

    end
end
disp(['Numero de erros de bit (Hard Decision): ', num2str(erroHard)]);
disp(['Numero de erros de bit (Soft Decision): ', num2str(erroSoft)]);

disp(['Numero de erros de bit (Hard Decision): ', num2str(numErrorsHard)]);
disp(['Taxa de erros de bit (Hard Decision): ', num2str(ratioHard)]);
disp(['Numero de erros de bit (Soft Decision): ', num2str(numErrorsSoft)]);
disp(['Taxa de erros de bit (Soft Decision): ', num2str(ratioSoft)]);
% disp(['Numero de erros de bit (Soft Decision): ', num2str(numErrorsUnCode)]);
% disp(['Taxa de erros de bit (Soft Decision): ', num2str(ratioUnCode)]);



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
title('Comparação de BER para Codificação LDPC');