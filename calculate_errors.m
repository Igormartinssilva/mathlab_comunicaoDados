%Cálculo de BER e FER
function [ber, fer] = calculate_errors(tx_bits, rx_bits, frame_size)
    errors = sum(tx_bits ~= rx_bits);
    ber = errors / length(tx_bits);
    fer = errors > 0; % FER é 1 se houve pelo menos 1 erro no quadro
end
