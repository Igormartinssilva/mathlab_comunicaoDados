%TRABALHO FINAL
%GRUPO 4 - IGOR MARTINS SILVA (00333069) E ARTUR RUIZ DE SOUZA (00334954)
close all;
clear;
clc;
%declara��o de vari�veis
mod_QPSK = 4; 
mod_QAM = 64;
n = 1296; %comprimento do c�digo
num_b = n; 
R = 2/3;
k = R*n; %bits de mensagem
bits_pari = n - k; %bits de paridade
frame_size = 2300 * 8;
Eb_N0_dB = -2:1:12; % Faixa de Eb/N0 em dB
Eb_N0_lin = 10 .^ (Eb_N0_dB / 10); % Faixa de Eb/N0 em linearizada
num_frames = 100; % N�mero de quadros simulados por Eb/N0

mensagem = randi(2,1,num_b)-1;
mensagemDemod = zeros(1, num_b);

%declara��o das matrizes e vetores
ber_qpsk = zeros(3, length(Eb_N0_lin));
fer_qpsk = zeros(3, length(Eb_N0_lin));
ber_qam = zeros(3, length(Eb_N0_lin));
fer_qam = zeros(3, length(Eb_N0_lin));


% QPSK ---------------------------------------------------------

%Vari�veis
Eb = sqrt(2); % Energia m�dia para quadratura
NP = Eb ./ (Eb_N0_lin); %vetor de pot�ncias do ru�do
NA = sqrt(NP); %vetor de amplitudes do ru�do
SI = zeros(1, num_b/2); 
SQ = zeros(1, num_b/2);

% Modula��o da mensagem em QPSK, construindo vetores para Q e I
cont = 1; % contador para evitar buracos nos vetores Q e I
for i = 1:2:num_b-1
    if (mensagem(i) == 1)
        SI(cont) = 1;
    else
        SI(cont) = -1;
    end
    if (mensagem(i+1) == 1)
        SQ(cont) = 1;
    else
        SQ(cont) = -1;
    end
    cont = cont + 1;
end

SMod = SI + 1i .* SQ; %Sinal resultante, modulado e complexo

for i = 1:length(Eb_N0_lin)
    N = NA(i)*complex(randn(1, num_b/2), randn(1, num_b/2))*sqrt(0.5); %vetor de ru�do complexo com desvio padr�o igual a uma posi��o do vetor NA
    r = SMod + N ; % vetor recebido
    cont = 1;
    for j = 1:2:num_b-1  %demodula��o do QPSK ap�s somar ru�do
        if (real(r(cont)) > 0) % valor do bit definido pelo quadrante em que o sinal Q+I se encontra
            mensagemDemod(j) = 1;
        else
            mensagemDemod(j) = 0;
        end
        if (imag(r(cont)) > 0)
            mensagemDemod(j+1) = 1;
        else
            mensagemDemod(j+1) = 0;
        end
        cont = cont + 1;
    end
    ber_qpsk(1, i) = sum(mensagem ~= mensagemDemod) / num_b; % contagem de erros e c�lculo do BER
end

%semilogy(Eb_N0_dB, ber_qpsk(1,:), 'r', 'LineWidth', 2, 'MarkerSize', 10);
%xlabel('Eb/N0 (dB)');
%ylabel('BER');
%legend('QPSK sem Cod');

%--------------------------------------LDPC encoder------------------------

% Criar matriz H na forma [P | I]
P = randi([0 1], bits_pari, k); % Parte aleat�ria (paridade)
I = eye(bits_pari);             % Matriz identidade
H = [P I];              % Combina��o de P e I

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
figure;
scatter(real(qam_symbols), imag(qam_symbols), 'filled');
xlabel('Eixo I (In-Phase)');
ylabel('Eixo Q (Quadrature)');
title('Constela��o 64-QAM');
grid on;
axis([-8 8 -8 8]); % Ajustar os limites para visualiza��o correta













