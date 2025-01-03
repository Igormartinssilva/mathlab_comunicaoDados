function H = qc_matrix_1296(n,r)
     % Gera uma matriz de paridade quase-c�clica esparsa v�lida
    % n: comprimento do c�digo
    % r: taxa do c�digo
    % Calcula o n�mero de linhas e colunas da matriz de paridade
    k = round(n * r);
    m = n - k; %bits de paridade

    % Define par�metros para a matriz quase-c�clica
    z = 54;  % Tamanho do bloco c�clico (deve ser divisor de m)
    nb = m / z; % N�mero de blocos de linha
    mb = n / z; % N�mero de blocos de coluna

    % Inicializa a matriz de paridade
    H = sparse(m, n);
 % Inicializa a matriz de deslocamento
 %int H_1296_2_3[8][24] = Matriz deslocamento para 1296 bits e razao 2/3
    shiftMatrix = [
        39, 31, 22, 43, -1, 40,  4, -1, 11, -1, -1, 50, -1, -1, -1,  6,  1,  0, -1, -1, -1, -1, -1, -1;
        25, 52, 41,  2,  6, -1, 14, -1, 34, -1, -1, -1, 24, -1, 37, -1, -1,  0,  0, -1, -1, -1, -1, -1;
        43, 31, 29,  0, 21, -1, 28, -1, -1,  2, -1, -1,  7, -1, 17, -1, -1, -1,  0,  0, -1, -1, -1, -1;
        20, 33, 48, -1,  4, 13, -1, 26, -1, -1, 22, -1, -1, 46, 42, -1, -1, -1, -1,  0,  0, -1, -1, -1;
        45,  7, 18, 51, 12, 25, -1, -1, -1, 50, -1, -1,  5, -1, -1, -1,  0, -1, -1, -1,  0,  0, -1, -1;
        35, 40, 32, 16,  5, -1, -1, 18, -1, -1, 43, 51, -1, 32, -1, -1, -1, -1, -1, -1, -1,  0,  0, -1;
         9, 24, 13, 22, 28, -1, -1, 37, -1, -1, 25, -1, -1, 52, -1, 13, -1, -1, -1, -1, -1, -1,  0,  0;
        32, 22,  4, 21, 16, -1, -1, -1, 27, 28, -1, 38, -1, -1, -1,  8,  1, -1, -1, -1, -1, -1, -1,  0;
    ];
    
    % Constr�i a matriz H a partir da matriz de deslocamento
    for i = 1:nb
        for j = 1:mb
            if shiftMatrix(i, j) >= 0
                H((i-1)*z+1:i*z, (j-1)*z+1:j*z) = circshift(eye(z), [0, shiftMatrix(i, j)]);
            end
        end
    end

    % Assegura que as �ltimas (N-K) colunas s�o invert�veis em GF(2)
    lastColumns = full(H(:, end-m+1:end)); % Converte para matriz densa
    while rank(lastColumns) < m
        % Regenera as �ltimas colunas se n�o forem de posto completo
        for i = 1:nb
            for j = mb-nb+1:mb
                if rand < 0.9   % Probabilidade de ter um bloco n�o-zero (ajuste conforme necess�rio)
                    block = circshift(eye(z), randi([0 z-1]));
                    lastColumns((i-1)*z+1:i*z, (j-1)*z+1:j*z) = block;
                end
            end
        end
    end
    H(:, end-m+1:end) = sparse(lastColumns); % Converte de volta para matriz esparsa
end