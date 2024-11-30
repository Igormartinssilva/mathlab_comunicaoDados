%Igor Martins Silva - 00333069
close all;
clear;

Nt = 4;%numero de antenas transmitindo
Nr = 4;%numero de antenas recebendo

tx_simb  = complex(2  * randi([0 1], Nt, 1) - 1); %vetor de simbolos baseados na quant. antenas
H = randn(Nr, Nt) + 1j * randn(Nr,Nt); %geraçao de matriz H aleatoria complexa

H_pseudo = pinv(H);% pseudo-inversa da H

detect_antes = H * tx_simb; %simbolos recebidos
detect_depois = H_pseudo * detect_antes; %simbolos dectados e achados pela tecnica Zero-Forcing


disp('Simbolos transmitidos:');
disp(tx_simb);

disp('Matriz H:');
disp(H);

disp('Pseudo-inversa de H:');
disp(H_pseudo);

disp('Simbolos recebidos (antes da detecção):');
disp(detect_antes);

disp('Simbolos recebidos (apos a detecção):');
disp(detect_depois);


