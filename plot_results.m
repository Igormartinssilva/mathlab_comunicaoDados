function plot_results(Eb_N0_dB, ber_qpsk, fer_qpsk, ber_qam, fer_qam)
    figure;
    semilogy(Eb_N0_dB, ber_qpsk(1, :), '-o', 'DisplayName', 'QPSK sem código');
    hold on;
    semilogy(Eb_N0_dB, ber_qpsk(2, :), '-x', 'DisplayName', 'QPSK com LDPC hard');
    semilogy(Eb_N0_dB, ber_qpsk(3, :), '-s', 'DisplayName', 'QPSK com LDPC soft');
    semilogy(Eb_N0_dB, ber_qam(1, :), '-^', 'DisplayName', '64-QAM sem código');
    semilogy(Eb_N0_dB, ber_qam(2, :), '-d', 'DisplayName', '64-QAM com LDPC hard');
    semilogy(Eb_N0_dB, ber_qam(3, :), '-p', 'DisplayName', '64-QAM com LDPC soft');
    xlabel('Eb/N0 (dB)');
    ylabel('BER');
    legend;
    grid on;
    title('Desempenho BER');
end
