function encoded_bits = ldpc_encode(message_bits, H, n, k, bits_pari)
    % Codificação LDPC sistemática
    % message_bits: vetor de k bits da mensagem
    % H: matriz de paridade (m x n)
    % n: comprimento total do código
    % k: número de bits da mensagem

    % Validar dimensões da matriz H
    if size(H, 1) ~= bits_pari || size(H, 2) ~= n
        error('As dimensões da matriz H não correspondem a m x n.');
    end

    % Garantir que H esteja na forma [P | I], logo separar as matrizes
    P = H(:, 1:k); % Parte que multiplica os bits de mensagem
    I = H(:, k+1:end); % Parte identidade para os bits de paridade

    if ~isequal(I, eye(bits_pari))
        error('A matriz H não está na forma [P | I].');
    end

    % Calcular os bits de paridade
    parity_bits = mod(P * message_bits, 2); % Multiplicação binária (P x mensagem)

    % Concatenar bits de mensagem e paridade
    encoded_bits = [message_bits; parity_bits];
end
