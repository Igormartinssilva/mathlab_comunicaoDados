function H = ldpcMatrizDeParidade(n, R)
  % Function to generate an LDPC parity-check matrix in systematic form
    % n: Length of the codeword
    % R: Rate of the code (e.g., 2/3)

    % Calculate the number of rows (m) and columns (n) of the parity-check matrix
    k = floor(n * R);  % Number of information bits
    m = n - k;         % Number of parity-check bits

    % Create an empty sparse logical matrix
    H = sparse(false(m, n));

    % Fill the matrix with random values while ensuring the last (N-K) columns
    % form an identity matrix
    for i = 1:m
        % Randomly select 'rowWeight' unique positions in the row to be ones
        % excluding the last (N-K) columns
        rowWeight = round(n / m);
        positions = randperm(k, rowWeight);
        H(i, positions) = true;
    end

    % Ensure the last (N-K) columns form an identity matrix
    H(:, k+1:end) = speye(m);
% Ensure each column has at least one '1'
    for j = 1:k
        if sum(H(:, j)) == 0
            % Randomly select a row and set the position to one
            rowIdx = randi(m);
            H(rowIdx, j) = true;
        end
    end
    % Convert to logical matrix
    H = logical(H);
end
