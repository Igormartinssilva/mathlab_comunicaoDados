function qpsk_symbols = manual_qpskmod(input_bits)
   % Modulação QPSK
    % input_bits: vetor de bits (codificado ou não)
    % Retorna: símbolos QPSK no plano complexo
    bits_per_symbol = 2; % Número de bits por símbolo
    num_b = length(input_bits);
    
    SI = zeros(1, num_b/2); % parte real do sinal modulado
    SQ = zeros(1, num_b/2); % parte imaginária do sinal modulado
    
    % Verificar se o número de bits é múltiplo de bits_per_symbol
    if mod(length(input_bits), bits_per_symbol) ~= 0
        error('O número de bits deve ser múltiplo de %d.', bits_per_symbol);
    end
    cont = 1; %contador para evitar "buracos" no vetor modulado
    for i = 1:2:num_b-1
        if (input_bits(i) == 1)
            SI(cont) = 1;
        else
            SI(cont) = -1;
        end
        if (input_bits(i+1) == 1)
            SQ(cont) = 1;
        else
            SQ(cont) = -1;
        end
        cont = cont + 1;
    end

    qpsk_symbols = SI + 1i .* SQ; 
end