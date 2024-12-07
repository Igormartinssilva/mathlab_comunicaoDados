function encoded_bits = ldpc_encode(message_bits, H, n, k, bits_pari)
    % Codifica��o LDPC sistem�tica
    % message_bits: vetor de k bits da mensagem
    % H: matriz de paridade (m x n)
    % n: comprimento total do c�digo
    % k: n�mero de bits da mensagem

    % Validar dimens�es da matriz H
    if size(H, 1) ~= bits_pari || size(H, 2) ~= n
        error('As dimens�es da matriz H n�o correspondem a m x n.');
    end

    % Garantir que H esteja na forma [P | I], logo separar as matrizes
    P = H(:, 1:k); % Parte que multiplica os bits de mensagem
    I = H(:, k+1:end); % Parte identidade para os bits de paridade

    if ~isequal(I, eye(bits_pari))
        error('A matriz H n�o est� na forma [P | I].');
    end

    % Calcular os bits de paridade
    parity_bits = mod(P * message_bits, 2); % Multiplica��o bin�ria (P x mensagem)

    % Concatenar bits de mensagem e paridade
    encoded_bits = [message_bits; parity_bits];
end
