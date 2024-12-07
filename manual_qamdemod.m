function demod_bits = manual_qamdemod(rx_symbols, M)
    % Demodulação manual para constelações QAM
    % rx_symbols: símbolos recebidos no plano complexo
    % M: ordem da modulação (ex.: 4 para QPSK, 16 para 16-QAM, etc.)
    
    % Verificar se M é potência de 2
    if mod(log2(M), 1) ~= 0
        error('M deve ser uma potência de 2');
    end

    % Criar constelação QAM
    sqrt_M = sqrt(M);
    I = repmat(-(sqrt_M-1):2:(sqrt_M-1), sqrt_M, 1);
    Q = repmat((sqrt_M-1):-2:-(sqrt_M-1), sqrt_M, 1)';
    constelation = I(:) + 1j * Q(:);
    
    % Normalizar a constelação para ter energia média unitária
    constelation = constelation / sqrt(mean(abs(constelation).^2));
    
    % Inicializar saída
    demod_bits = zeros(length(rx_symbols), log2(M));
    
    % Para cada símbolo recebido, encontrar o mais próximo na constelação
    for k = 1:length(rx_symbols)
        [~, idx] = min(abs(rx_symbols(k) - constelation)); % Índice do símbolo mais próximo
        demod_bits(k, :) = de2bi(idx-1, log2(M), 'left-msb'); % Converte índice para bits
    end
    
    % Converter para vetor de bits
    demod_bits = demod_bits(:);
end
