function qam_symbols = manual_qammod(input_bits, mod_QAM)
   % Modulação 64-QAM
    % input_bits: vetor de bits (codificado ou não)
    % Retorna: símbolos QAM no plano complexo
    bits_per_symbol = log2(mod_QAM); % Número de bits por símbolo (6 bits para 64-QAM)

    % Verificar se o número de bits é múltiplo de bits_per_symbol
    if mod(length(input_bits), bits_per_symbol) ~= 0
        error('O número de bits deve ser múltiplo de %d.', bits_per_symbol);
    end

    % Número de símbolos
    num_symbols = length(input_bits) / bits_per_symbol;

    % Agrupar bits em blocos de 6
    bit_groups = reshape(input_bits, bits_per_symbol, num_symbols)';

    % Converter grupos de bits para índices decimais
    decimal_indices = manual_bi2de(bit_groups, true); % Usando função manual_bi2de

    % Gerar constelação 64-QAM
    sqrt_M = sqrt(mod_QAM); % Tamanho da grade
    I = repmat(-(sqrt_M-1):2:(sqrt_M-1), sqrt_M, 1); % Eixo I (in-phase)
    Q = repmat((sqrt_M-1):-2:-(sqrt_M-1), sqrt_M, 1)'; % Eixo Q (quadrature)
    constelation = I(:) + 1j * Q(:); % Constelação no plano complexo
    
    qam_symbols = constelation;
    
    % Normalizar a constelação (energia média unitária)
     %constelation = constelation / sqrt(mean(abs(constelation).^2));

    % Mapear os índices decimais para símbolos QAM
    %qam_symbols = constelation(decimal_indices + 1);
    
end
