%Igor Martins Silva - 00333069
close all;
clear;

Nt = 4;
Ns = 4;

tx_simb  = 2  * randi([0 1], Nt, 1) - 1
H = randn(NR, NT) + 1j * randn(NR,NT)

H_pseudo = pinv(H);

detect_antes =  tx_simb * H_pseudo;

detect = sing(real(detect_antes));

disp('Símbolos transmitidos:');
disp(tx_simb);

disp('Matriz H:');
disp(H);

disp('Pseudo-inversa de H:');
disp(H_pseudo);

disp('Símbolos recebidos (antes da detecção):');
disp(detect_antes);

disp('Símbolos detectados (após a detecção):');
disp(detect);
