% Simulation parameters
N = 1e6;
snrRange = 0:2:30;
m = 1; % Number of samples for each bit (for OOK and PRK)

% Generate random binary data vector
rand = randi([0, 1], 1, N);

% Modulation: ASK (OOK)
askModulated = rand;

%FSK special case: Orthogonal-FSK modulation
ones=find(rand);     %find all elements that are equal to one in the rand sequence
zeros=find(~rand);   %find all elements that are equal to zero in the rand sequence
fskModulated(ones)=i;  %replace all the ones with i which is an imaginary number
fskModulated(zeros)=1; %replace all the zeros with 1

% Modulation: PSK (PRK)
pskModulated = (rand.*2) - 1;

% Initialize BER matrices
ber_ask = zeros(1, length(snrRange));
ber_fsk = zeros(1, length(snrRange));
ber_psk = zeros(1, length(snrRange));
berAsk2 = zeros(1, length(snrRange));
berFsk2 = zeros(1, length(snrRange));
berPsk2 = zeros(1, length(snrRange));
berQam = zeros(1, length(snrRange));
% modulating using built in function 
%random_data=randi([0 16-1],1,N);              %Generate random data vector(1 x n)
ask2=genqammod(rand,[0  1]);   %OOK modulation using built in functions  
psk2=pskmod(rand,2);           %PRK modulation using built in functions  
fsk2=genqammod(rand,[1  1i]);  %FSK modulation using built in functions 
qam=qammod(rand,16);           %QAM modulation using built in functions    


% Perform simulations for each SNR
for snrIndex = 1:length(snrRange)
    snr = snrRange(snrIndex);
    
    % Apply noise to modulated signals
    % Calculate noise power based on SNR
    noise = sqrt( 1 / (10^(snr/10))) * randn(1, N);

    askRecieved2 = ask2 + noise;
    fskRecieved2 = fsk2 + noise;
    pskRecieved2 = psk2 + noise;
    qamRecevied = qam + noise;
    
    askReceived = askModulated + noise;
    fskReceived = fskModulated + noise;
    pskReceived = pskModulated + noise;
    
    % Demodulation using built-in functions
    askDetected2 = genqamdemod(askRecieved2, [0 1]);
    fskDetected2 = genqamdemod(fskRecieved2, [1 1i]);
    pskDetected2 = pskdemod(pskRecieved2, 2);
    qamDetected2 = qamdemod(qamRecevied,16);

    % Decide whether the received signal is '1' or '0' for each modulation scheme
    askDetected = real(askReceived) >= 0.5;
    fskDetected = real(fskReceived) < imag(fskReceived);
    pskDetected = real(pskReceived) >= 0;
    
    % Calculate number of bit errors
    numErrors_ask = biterr(rand, askDetected);
    numErrors_fsk = biterr(rand, fskDetected);
    numErrors_psk = biterr(rand, pskDetected);
    
    % Calculate bit error rate (BER)
    ber_ask(snrIndex) = numErrors_ask / N;
    ber_fsk(snrIndex) = numErrors_fsk / N;
    ber_psk(snrIndex) = numErrors_psk / N;

    % Calculate number of bit errors for built-in
    numErrorsAsk2 = biterr(rand, askDetected2);
    numErrorsFsk2 = biterr(rand, fskDetected2);
    numErrorsPsk2 = biterr(rand, pskDetected2);
    numErrorsQam = biterr(rand, qamDetected2);
    
    % Calculate bit error rate (BER) for built-in
    berAsk2(snrIndex) = numErrorsAsk2 / N;
    berFsk2(snrIndex) = numErrorsFsk2 / N;
    berPsk2(snrIndex) = numErrorsPsk2 / N;
    berQam(snrIndex) = numErrorsQam / N;
end
% Find SNR value for nearly error-free transmission
snr_ask_min = snrRange(find(ber_ask <= 1e-5, 1));
snr_fsk_min = snrRange(find(ber_fsk <= 1e-5, 1));
snr_psk_min = snrRange(find(ber_psk <= 1e-5, 1));
snrAskMin2 = snrRange(find(berAsk2 <= 1e-5, 1));
snrFskMin2 = snrRange(find(berFsk2 <= 1e-5, 1));
snrPskMin2 = snrRange(find(berPsk2 <= 1e-5, 1));
snrQAMMin = snrRange(find(berQam <= 1e-5, 1));

disp("SNR for nearly error-free transmission:");
disp("ASK (OOK): " + snr_ask_min + " dB");
disp("FSK (Orthogonal-FSK): " + snr_fsk_min + " dB");
disp("PSK (PRK): " + snr_psk_min + " dB");
disp("ASK (OOK) Built-in: " + snrAskMin2 + " dB");
disp("FSK (Orthogonal-FSK) Built-in: " + snrFskMin2 + " dB");
disp("PSK (PRK) Built-in: " + snrPskMin2 + " dB");
disp("QAM Built-in: " + snrQAMMin + " dB");

% Plot the BER curves for different modulation schemes
figure;
semilogy(snrRange, ber_ask, '-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(snrRange, ber_fsk, '-s', 'LineWidth', 2, 'MarkerSize', 8);
semilogy(snrRange, ber_psk, '-d', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(snrRange, berAsk2, '--+', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(snrRange, berFsk2, '--*', 'LineWidth', 2, 'MarkerSize', 8);
semilogy(snrRange, berPsk2, '--x', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(snrRange, berQam, '--', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
title('Bit Error Rate vs SNR');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
legend('ASK (OOK)', 'FSK (Orthogonal-FSK)', 'PSK (PRK)','ASK (OOK) Built-in', 'FSK (Orthogonal-FSK) Built-in', 'PSK (PRK) Built-in','QAM Built-in');