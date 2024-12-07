function decimal = manual_bi2de(binary_matrix, msb_first)
    % Convertendo matriz binária para decimal
    if nargin < 2
        msb_first = true;
    end

    % Inverter os bits se LSB primeiro
    if ~msb_first
        binary_matrix = fliplr(binary_matrix);
    end

    % Converter para decimal
    [num_rows, num_bits] = size(binary_matrix);
    decimal = zeros(num_rows, 1);
    for i = 1:num_bits
        decimal = decimal + binary_matrix(:, i) * 2^(num_bits - i);
    end
end
