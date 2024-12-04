%TRABALHO FINAL
%GRUPO 4 - IGOR MARTINS SILVA (00333069) E ARTUR RUIZ DE SOUZA (00334954)
close all;
clear;

%declaração de variáveis
mod_QPSK = 4; 
mod_QAM = 64;
n = 1296;
num_b = n;
R = 2/3;
frame_size = 2300 * 8;
Eb_N0_dB = -2:1:12; % Faixa de Eb/N0 em dB
Eb_N0_lin = 10 .^ (Eb_N0_dB / 10); % Faixa de Eb/N0 em linearizada
num_frames = 100; % Número de quadros simulados por Eb/N0
mensagem = randi(2,1,num_b)-1;
mensagemDemod = zeros(1, num_b);

%declaração das matrizes e vetores
ber_qpsk = zeros(3, length(Eb_N0_lin));
fer_qpsk = zeros(3, length(Eb_N0_lin));
ber_qam = zeros(3, length(Eb_N0_lin));
fer_qam = zeros(3, length(Eb_N0_lin));



% QPSK ---------------------------------------------------------

%Variáveis
Eb = sqrt(2); % Energia média para quadratura
NP = Eb ./ (Eb_N0_lin); %vetor de potências do ruído
NA = sqrt(NP); %vetor de amplitudes do ruído
SI = zeros(1, num_b/2); 
SQ = zeros(1, num_b/2);

% Modulação da mensagem em QPSK, construindo vetores para Q e I
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
    N = NA(i)*complex(randn(1, num_b/2), randn(1, num_b/2))*sqrt(0.5); %vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    r = SMod + N ; % vetor recebido
    cont = 1;
    for j = 1:2:num_b-1  %demodulação doQPSK após somar ruído
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
    ber_qpsk(1, i) = sum(mensagem ~= mensagemDemod) / num_b; % contagem de erros e cálculo do BER
end

semilogy(Eb_N0_dB, ber_qpsk(1,:), 'r', 'LineWidth', 2, 'MarkerSize', 10);
xlabel('Eb/N0 (dB)');
ylabel('BER');
legend('QPSK sem Cod');
























