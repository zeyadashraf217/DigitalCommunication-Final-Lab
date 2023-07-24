clear all;
clc;
% Step 1: Simulation parameters
numBits = 1e5;
snrRange = 0:2:30;
m = [10, 20, 100]; % Number of samples for s1(t)
samplingInstant = 20; % Sampling instant
s1Amplitude = 1;
s2Amplitude = 0;

% Step 2: Generate random binary data vector
bits = randi([0, 1], 1, numBits);

% Step 3: Represent each bit with proper waveform
s1 = cell(length(m), 1);
for i = 1:length(m)
    s1{i} = s1Amplitude * ones(1, m(i));
end
s2 = s2Amplitude * ones(1, max(m)); % Use the maximum number of samples for s2
waveform = kron(bits, s1{2});
%waveform = [];
%for i = 1:numBits
%    if bits(i) == 0
%        waveform = [waveform, s2];
%    else
%        waveform = [waveform, s1{2}];
%    end
% end


% Step 4: Apply noise to samples
ber_matched = zeros(length(m), length(snrRange));
ber_correlator = zeros(length(m), length(snrRange));

for mIndex = 1:length(m)
    for snrIndex = 1:length(snrRange)
        
        snr = snrRange(snrIndex);
        % Add white Gaussian noise to the waveform
        rxSequence = awgn(waveform, snr, 'measured');

        % Step 5: Apply convolution process in the receiver (Matched Filter)
        s1_minus_s2 = s1{mIndex}(1:length(s1{mIndex})) - s2(1:length(s1{mIndex}));
        MatchedOutput = conv(rxSequence, fliplr(s1_minus_s2));
        matchedSamples = MatchedOutput(m(mIndex):samplingInstant:end);
        

        % Step 6: Decide whether the Rx_sequence is '1' or '0' by comparing with threshold (Matched Filter)
        threshold_matched = 0.5 * max(matchedSamples);
        detectedBits_matched = matchedSamples > threshold_matched;

        % Step 7: Compare the original bits with the detected bits and calculate number of errors (Matched Filter)
        numErrors_matched = biterr(bits, detectedBits_matched);

        % Step 8: Save the probability of error of each SNR in matrix, BER (Matched Filter)
        ber_matched(mIndex, snrIndex) = numErrors_matched / numBits;

        % Step 9: Apply correlation process in the receiver (Correlator)
        noise_power = 1/snr;
        noise = sqrt(noise_power) * randn(1, numBits);
        correlatorOutput=bits .* noise;
        threshold_Corr = sum(correlatorOutput) / length(correlatorOutput);
        correlatorSamples = correlatorOutput(m(mIndex):samplingInstant:end);
        % Step 10: Decide whether the Rx_sequence is '1' or '0' by comparing with threshold (Correlator)
        detectedBits_correlator = correlatorSamples > threshold_Corr;

        % Step 11: Compare the original bits with the detected bits and calculate number of errors (Correlator)
        numErrors_correlator = biterr( (bits(1:length(detectedBits_correlator))) , detectedBits_correlator);
        
        % Step 12: Save the probability of error of each SNR in matrix, BER (Correlator)
        ber_correlator(mIndex, snrIndex) = numErrors_correlator / numBits;
    end
end

% Step 13: Plot BER vs SNR (Matched Filter and Correlator)
figure;
for mIndex = 1:length(m)
    semilogy(snrRange, ber_matched(mIndex, :), '-o', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
end
for mIndex = 1:length(m)
    semilogy(snrRange, ber_correlator(mIndex, :), '--o', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
end
grid on;
title('Bit Error Rate vs SNR');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
legend('Matched Filter (m=10)', 'Matched Filter (m=20)', 'Matched Filter (m=100)', 'Correlator (m=10)', 'Correlator (m=20)', 'Correlator (m=100)');

% Step 14: Calculation of transmitted signal power
transmittedPower = mean(waveform .^ 2);
fprintf('Transmitted Signal Power: %.2f\n', transmittedPower);

% Step 15: At which value of SNR the system is nearly without error (for the given frame)?
thresholdBER = 1e-5;
snrWithoutError = zeros(1, length(m));
for mIndex = 1:length(m)
    [~, idx] = min(ber_matched(mIndex, :));
    snrWithoutError(mIndex) = snrRange(idx);
end
fprintf('SNR at which the system is nearly without error (BER <= %.0e):\n', thresholdBER);
fprintf('m=10: SNR = %d dB\n', snrWithoutError(1));
fprintf('m=20: SNR = %d dB\n', snrWithoutError(2));
fprintf('m=100: SNR = %d dB\n', snrWithoutError(3));

% Step 16: Additional part to calculate BER without matched filter and correlator
BER = zeros(1, length(snrRange));
iterations = 1;

for temp = 1:length(snrRange)
    SNR = 10^(snrRange(temp)/10);
    % Average power of transmitted signal / SNR
    noise_power = 1/SNR;
    noise = sqrt(noise_power) * randn(1, numBits);
    num_errors = 0;

    for i = 1:iterations
        data = randi([0 1], 1, numBits);
        % Apply noise to signal
        received_signal = data + noise;
        % Decide whether the received signal is '1' or '0'
        detected_data = (received_signal >= 0.5);
        errors = biterr(data, detected_data);
        num_errors = num_errors + errors;
    end
    % Calculate bit error rate (BER)
    BER(temp) = num_errors / (iterations*numBits);
end

% Step 17: Plot the additional BER curve
figure;
semilogy(snrRange, BER, '-o', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
title('Bit Error Rate vs SNR (Without Matched Filter and Correlator)');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
legend('Without Matched Filter and Correlator');
