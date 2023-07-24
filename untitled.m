% Number of bits/SNR
n=1e3;
SNR=0:2:30;
snr=10.^(SNR/10);
%initializing empty matrices to be filled with the sequences
OOK=zeros(1,16);
PRK=zeros(1,16);
ASK=zeros(1,16);
PSK=zeros(1,16);
FSK=zeros(1,16);
OOK_mat=zeros(1,16);
PRK_mat=zeros(1,16);
ASK_mat=zeros(1,16);
PSK_mat=zeros(1,16);
FSK_mat=zeros(1,16);
QAM_mat=zeros(1,16);
% Transmitter side variables
%Generate random binary data vector(1 x n)
data=randi([0,1],1,n);
%we will work on ASK PSK FSK modulation according to the specified cases inthe manual
%ASK special case: OOK_Modulation(0,1)
signal_OOK=data; %No change in the bits will be required
%PSK special case: PRK modulation(-1,1)
signal_PRK1=(data.*2)-1; %represent the 1 by 1 and the 0 bit by -1 using providedformula
signal_PRK=complex(signal_PRK1,0);
%ASK
signal_ASK = data;
signal_ASK(data == 0) = -1;
signal_ASK(data == 1) = 1;
%PSK
signal_PSK=data;
signal_PSK(data== 1) = 1i;
signal_PSK(data == 0) = -1i;
%FSK special case: Orthogonal-FSK modulation
signal_FSK = data;
% ones=find(data); %find all elements that are equal to one in the bits sequence
% zeros=find(~data); %find all elements that are equal to zero in the bitssequence
% signal_FSK(ones)=i; %replace all the ones with i which is an imaginary number
% signal_FSK(zeros)=1; %replace all the zeros with 1
for i = 1:n
 if signal_FSK(i) == 0
 signal_FSK(i) = 1;
 else
 signal_FSK(i) = 1j;
 end
end
% modulating using matlab built in function
OOK_mat_mod=genqammod(data,[0 1]); %OOK modulation using built infunctions
PRK_mat_mod=pskmod(data,2); %PRK modulation using built in functions, pulse width = 2 (1,-1)
ASK_mat_mod=pammod(data,2); %ASK modulation using built infunctions
PSK_mat_mod=pskmod(data,2); %PSK modulation using built in functions
FSK_mat_mod=genqammod(data,[1 1i]); %FSK modulation using built infunctions
rbdv=randi([0 15],1,n); %Generate random data vector(1 x n)
QAM_mat_mod=qammod(rbdv,16); %QAM modulation using built infunctions , modulation order = 16
%Applying noise to the bit stream with different values of SNR,receiving
for i=1:16
 %applying the noise with mean = 0 and varience =1
 AvgPower_OOK=mean(signal_OOK.^2);
 %noise=(sqrt(AvgPower_OOK/(2*snr(i))))*(randn(1,n)+(1j*randn(1,n)));
 OOK_sequence=awgn(signal_OOK,SNR(i),'measured');

 AvgPower_PRK=mean(signal_PRK.^2);
 %noise=(sqrt(AvgPower_PRK/(2*snr(i))))*(randn(1,n)+(1j*randn(1,n)));
 PRK_sequence=awgn(signal_PRK,SNR(i),'measured');

 AvgPower_ASK=mean(signal_ASK.^2);
 %noise=(sqrt(AvgPower_ASK/(2*snr(i))))*(randn(1,n)+(1j*randn(1,n)));
 ASK_sequence=awgn(signal_ASK,SNR(i),'measured');

 AvgPower_PSK=mean(signal_PSK.^2);
 %noise=(sqrt(AvgPower_PSK/(2*snr(i))))*(randn(1,n)+(1j*randn(1,n)));
 PSK_sequence=awgn(signal_PSK,SNR(i),'measured');

 AvgPower_FSK=mean(signal_FSK.^2);
 %noise=(sqrt(AvgPower_FSK /(2*snr(i))))*(randn(1,n)+(1j*randn(1,n)));
 FSK_sequence=awgn(signal_FSK,SNR(i),'measured');


 %aplying noise in built-in functions

%noise=(sqrt(mean(abs(OOK_mat_mod).^2)/(2*snr(i))))*(randn(1,n)+(1j*randn(1,n)));
 OOK_sequence_mat=awgn(OOK_mat_mod,SNR(i),'measured');


%noise=(sqrt(mean(abs(PRK_mat_mod).^2)/(2*snr(i))))*(randn(1,n)+(1j*randn(1,n)));
 PRK_sequence_mat=awgn(PRK_mat_mod,SNR(i),'measured');


%noise=(sqrt(mean(abs(ASK_mat_mod).^2)/(2*snr(i))))*(randn(1,n)+(1j*randn(1,n)));
 ASK_sequence_mat=awgn(ASK_mat_mod,SNR(i),'measured');


%noise=(sqrt(mean(abs(PSK_mat_mod).^2)/(2*snr(i))))*(randn(1,n)+(1j*randn(1,n)));
 PSK_sequence_mat=awgn(PSK_mat_mod,SNR(i),'measured');


%noise=(sqrt(mean(abs(FSK_mat_mod).^2)/(2*snr(i))))*(randn(1,n)+(1j*randn(1,n)));
 noise = awgn(FSK_mat_mod,SNR(i),'measured');
 FSK_sequence_mat=FSK_mat_mod+noise;


%noise=(sqrt(mean(abs(QAM_mat_mod).^2)/(2*snr(i))))*(randn(1,n)+(1j*randn(1,n)));
 QAM_sequence_mat=awgn(QAM_mat_mod,SNR(i),'measured');

 %Decide whether the Rx_sequence is ‘1’ or ‘0’ by comparing each bit with athreshold
 OOK_demod=(real(OOK_sequence) >=0.5);%0.5 for OOK
 PRK_demod=(real(PRK_sequence) >=0);%0 for PRK
 ASK_demod=(real(ASK_sequence) >=0);%0 for ASK
 PSK_demod=(imag(PSK_sequence) >=0);%0 for PSK
 FSK_demod=(real(FSK_sequence)<imag(FSK_sequence));%in FSK if real >imaginary component then it equal zero else equals one

 %demodulation using built in functions
 OOK_mat_demod = genqamdemod(OOK_sequence_mat,[0 1]);
 PRK_mat_demod = pskdemod(PRK_sequence_mat,2);
 ASK_mat_demod = pamdemod(ASK_sequence_mat,2);
 PSK_mat_demod = pskdemod(PSK_sequence_mat,2);
 FSK_mat_demod = genqamdemod(FSK_sequence_mat,[1 1i]);
 rbdv_mat_demod = qamdemod(QAM_sequence_mat,16);

 %Calculate #of errors (BER) %[NUMBER,RATIO] = BITERR(X,Y)
 [~,OOK(i)]=biterr(data,OOK_demod);
 [~,PRK(i)]=biterr(data,PRK_demod);
 [~,ASK(i)]=biterr(data,ASK_demod);
 [~,PSK(i)]=biterr(data,PSK_demod);
 [~,FSK(i)]=biterr(data,FSK_demod);

 %for modulated using matlab
 [~,OOK_mat(i)] = biterr(data,OOK_mat_demod);
 [~,PRK_mat(i)] = biterr(data,PRK_mat_demod);
 [~,ASK_mat(i)] = biterr(data,ASK_mat_demod);
 [~,PSK_mat(i)] = biterr(data,PSK_mat_demod);
 [~,FSK_mat(i)] = biterr(data,FSK_mat_demod);
 [~,QAM_mat(i)] = symerr(rbdv,rbdv_mat_demod); %calculate symbol error rate
end
% Ploting curves
figure
semilogy(SNR,OOK,'r-*','LineWidth',2)
hold on;
semilogy(SNR,OOK_mat,'m--o','LineWidth',2)
title('BER vs. SNR (OOK) ')
ylabel('BER')
xlabel('SNR')
legend('OOK','OOK-Built-in function')
grid on;
figure
semilogy(SNR,PRK,'g-*','LineWidth',2)
hold on;
semilogy(SNR,PRK_mat,'k--+','LineWidth',2)
title('BER vs. SNR (PRK) ')
ylabel('BER')
xlabel('SNR')
legend('PRK','PRK-Built-in function')
grid on;
figure
semilogy(SNR,FSK,'b-*','LineWidth',2)
hold on;
semilogy(SNR,FSK_mat,'y--s','LineWidth',2)
title('BER vs. SNR (FSK) ')
ylabel('BER')
xlabel('SNR')
legend('FSK','FSK-Built-in function')
grid on;
figure
semilogy(SNR,ASK,'b-*','LineWidth',2)
hold on;
semilogy(SNR,ASK_mat,'y--s','LineWidth',2)
title('BER vs. SNR (ASK) ')
ylabel('BER')
xlabel('SNR')
legend('ASK','ASK-Built-in function')
grid on;
figure
semilogy(SNR,PSK,'b-*','LineWidth',2)
hold on;
semilogy(SNR,PSK_mat,'y--s','LineWidth',2)
title('BER vs. SNR (PSK) ')
ylabel('BER')
xlabel('SNR')
legend('PSK','PSK-Built-in function')
grid on;
figure
semilogy(SNR,QAM_mat,'c-p','LineWidth',2)
title('BER vs. SNR (QAM)')
ylabel('BER')
xlabel('SNR')
grid on;
figure
semilogy(SNR,OOK,'m--o','LineWidth',2)
hold on;
semilogy(SNR,PRK,'k--+','LineWidth',2)
hold on;
semilogy(SNR,FSK,'y--s','LineWidth',2)
title('BER vs. SNR (All Manual)')
ylabel('BER')
xlabel('SNR')
legend('OOK','PRK','FSK')
grid on;