function demod_bits = manual_qpskdemod(rx_symbols)
    % Demodulação manual para constelações QPSK
    % rx_symbols: símbolos recebidos no plano complexo
    % Retorna: símbolos QPSK no plano complexo

    cont = 1; %contador para evitar "buracos" no vetor demodulado
    num_b = length(rx_symbols)*2;
    mensagemDemod = zeros(1, num_b);
    for j = 1:2:num_b-1 
        if (real(rx_symbols(cont)) > 0) % valor do bit definido pelo quadrante em que o sinal Q+I se encontra
            mensagemDemod(j) = 1;
        else
            mensagemDemod(j) = 0;
        end
        if (imag(rx_symbols(cont)) > 0)
            mensagemDemod(j+1) = 1;
        else
            mensagemDemod(j+1) = 0;
        end
        cont = cont + 1;
    end
    demod_bits = mensagemDemod;
end