function demod_bits = manual_qamdemod(rx_symbols, M)
    % Demodula��o manual para constela��es QAM
    % rx_symbols: s�mbolos recebidos no plano complexo
    % M: ordem da modula��o (ex.: 4 para QPSK, 16 para 16-QAM, etc.)
    
    % Verificar se M � pot�ncia de 2
    if mod(log2(M), 1) ~= 0
        error('M deve ser uma pot�ncia de 2');
    end

    % Criar constela��o QAM
    sqrt_M = sqrt(M);
    I = repmat(-(sqrt_M-1):2:(sqrt_M-1), sqrt_M, 1);
    Q = repmat((sqrt_M-1):-2:-(sqrt_M-1), sqrt_M, 1)';
    constelation = I(:) + 1j * Q(:);
    
    % Normalizar a constela��o para ter energia m�dia unit�ria
    constelation = constelation / sqrt(mean(abs(constelation).^2));
    
    % Inicializar sa�da
    demod_bits = zeros(length(rx_symbols), log2(M));
    
    % Para cada s�mbolo recebido, encontrar o mais pr�ximo na constela��o
    for k = 1:length(rx_symbols)
        [~, idx] = min(abs(rx_symbols(k) - constelation)); % �ndice do s�mbolo mais pr�ximo
        demod_bits(k, :) = de2bi(idx-1, log2(M), 'left-msb'); % Converte �ndice para bits
    end
    
    % Converter para vetor de bits
    demod_bits = demod_bits(:);
end
