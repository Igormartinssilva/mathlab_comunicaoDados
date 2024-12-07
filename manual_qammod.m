function qam_symbols = manual_qammod(input_bits, mod_QAM)
   % Modula��o 64-QAM
    % input_bits: vetor de bits (codificado ou n�o)
    % Retorna: s�mbolos QAM no plano complexo
    bits_per_symbol = log2(mod_QAM); % N�mero de bits por s�mbolo (6 bits para 64-QAM)

    % Verificar se o n�mero de bits � m�ltiplo de bits_per_symbol
    if mod(length(input_bits), bits_per_symbol) ~= 0
        error('O n�mero de bits deve ser m�ltiplo de %d.', bits_per_symbol);
    end

    % N�mero de s�mbolos
    num_symbols = length(input_bits) / bits_per_symbol;

    % Agrupar bits em blocos de 6
    bit_groups = reshape(input_bits, bits_per_symbol, num_symbols)';

    % Converter grupos de bits para �ndices decimais
    decimal_indices = manual_bi2de(bit_groups, true); % Usando fun��o manual_bi2de

    % Gerar constela��o 64-QAM
    sqrt_M = sqrt(mod_QAM); % Tamanho da grade
    I = repmat(-(sqrt_M-1):2:(sqrt_M-1), sqrt_M, 1); % Eixo I (in-phase)
    Q = repmat((sqrt_M-1):-2:-(sqrt_M-1), sqrt_M, 1)'; % Eixo Q (quadrature)
    constelation = I(:) + 1j * Q(:); % Constela��o no plano complexo
    
    qam_symbols = constelation;
    
    % Normalizar a constela��o (energia m�dia unit�ria)
     %constelation = constelation / sqrt(mean(abs(constelation).^2));

    % Mapear os �ndices decimais para s�mbolos QAM
    %qam_symbols = constelation(decimal_indices + 1);
    
end
