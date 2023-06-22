%------------Designing an AM radio---------%

subplot(4,3,1);
[S1,fs1]= audioread("StarWars3.wav");
% sound(S1,fs1);
S1=resample(S1,3,1); 
fs1=fs1*3;
t1=1:3*fs1;
plot(t1/fs1,S1);
title("Plot of StarWars3.wav");
xlabel 'Time(sec)';
ylabel 'S1(t)';
grid on;

subplot(4,3,2);
fm1=3000;
S1=lowpass(S1,fm1,fs1);
S1_f=fftshift(fft(S1));
n1=length(S1_f);
freq=(-fs1/2):(fs1/(n1-1)):(fs1/2);
plot(freq,(abs(S1_f)/n1));
title("Frequency response of StarWars3.wav");
xlabel 'Frequency(Hz)';
ylabel 'S1(f)';
grid on;

subplot(4,3,5);
[S2,fs2]= audioread("CantinaBand3.wav");
% sound(S2,fs2);
S2=resample(S2,3,1);
fs2=fs2*3;
t2=1:3*fs2;
plot(t2/fs2,S2);
title("Plot of CantinaBand3.wav");
xlabel 'Time(sec)';
ylabel 'S2(t)';
grid on;

subplot(4,3,6);
fm2=3000;
S2=lowpass(S2,fm2,fs2);
S2_f=fftshift(fft(S2));
n2=length(S2_f);
freq1=(-fs2/2):(fs2/(n2-1)):(fs2/2);
plot(freq1,(abs(S2_f)/n2));
title("Frequency response of CantinaBand3.wav");
xlabel 'Frequency(Hz)';
ylabel 'S2(f)';
grid on;

fprintf("For StarWars3.wav\n");
disp("Sampling frequency: "+fs1+"Hz");
disp("Bandwidth: "+obw(S1,fs1)+"Hz");

fprintf("For CantinaBand3.wav\n");
disp("Sampling frequency: "+fs2+"Hz");
disp("Bandwidth: "+obw(S2,fs2)+"Hz");

fc1=10000;
fc2=25000;

A1=0;
A2=0;

%-----------------Transmitter-------------------%
S3=transmitter(S1,S2,fc1,fc2,t1,t2,fs1,fs2,A1,A2); 
subplot(4,3,9);
plot(t1/fs1,S3); %-----
title("Transmitted signal(S3)");
xlabel 'Time(sec)';
ylabel 'S3(t)';
grid on;    

subplot(4,3,10);
S3_f=fftshift(fft(S3));
n4=length(S3_f);
freq3=(-fs1/2):(fs1/(n4-1)):(fs1/2);
plot(freq3,(abs(S3_f)/n4));
title("Frequency response of S3");
xlabel 'Frequency(Hz)';
ylabel 'S3(F)';
grid on; 

%-----------------Channel------------------------%
[S3,snr]=channel(S3);   

%----------------Reciever------------------------%
prompt = "Select channel number ";
ch_number = input(prompt);
recieved_signal=reciever(ch_number,S3,fm1,fm2,fs1,fs2,t1,t2,fc1,fc2,A1,A2,snr); 

% --------------------------Channel---------------------%
function [add_noise,snr]=channel(am_signal)
    snr=10;
    add_noise=awgn(am_signal,snr,'measured');
    
end

% --------------------------Reciever---------------------%
function recieve_signal = reciever(channel_number,am_signal,fm1,fm2,fs1,fs2,t1,t2,fc1,fc2,A1,A2,noise)
    
    if(channel_number==1)
        subplot(4,3,11);
        [recieve_signal1,recieve_signal]=demod(am_signal,fc1,fm1,fs1,t1,A1);
        plot(t1/fs1,recieve_signal);
        title("Demodulated Signal-1");
        xlabel 'Time(sec)';
        ylabel 'Amplitude';
        grid on; 
        sound(recieve_signal,fs2);

        subplot(4,3,12);
        recieve_signal_f=fftshift(fft(recieve_signal));
        n5=length(recieve_signal_f);
        freq4=(-fs1/2):(fs1/(n5-1)):(fs1/2);
        plot(freq4,(abs(recieve_signal_f)/n5));
        title("Frequency response of Demodulated (S4)");
        xlabel 'Frequency(Hz)';
        ylabel 'S4(F)';
        grid on; 


    else
        subplot(4,3,11);
        [recieve_signal1,recieve_signal]=demod(am_signal,fc2,fm2,fs2,t2,A2);
        plot(t2/fs2,recieve_signal);
        title("Demodulated Signal-2");
        xlabel 'Time(sec)';
        ylabel 'Amplitude';
        grid on; 
        sound(recieve_signal,fs2);

        subplot(4,3,12);
        recieve_signal_f=fftshift(fft(recieve_signal));
        n5=length(recieve_signal_f);
        freq4=(-fs2/2):(fs2/(n5-1)):(fs2/2);
        plot(freq4,(abs(recieve_signal_f)/n5));
        title("Frequency response of Demodulated (S4)");
        xlabel 'Frequency(Hz)';
        ylabel 'S4(F)';
        grid on; 
    end
end

% --------------------------Transmitter---------------------%
function resultant_signal = transmitter(S1,S2,fc1,fc2,t1,t2,fs1,fs2,A1,A2)

    S1_c=cos(2*pi*(fc1/fs1)*t1);
    S1_modulated1=modul(S1,S1_c,A1);
    
    subplot(4,3,3);
    plot(t1/fs1,S1_modulated1);
    title("Modulated S1(t)");
    xlabel 'Time(sec)';
    ylabel 'Amplitude';
    grid on;

    subplot(4,3,4);
    S1_modulated1_f=fftshift(fft(S1_modulated1));
    n2=length(S1_modulated1_f);
    freq1=(-fs1/2):(fs1/(n2-1)):(fs1/2);
    plot(freq1,(abs(S1_modulated1_f)/n2));
    title("Frequency response of modulated S1(t)");
    xlabel 'Frequency(Hz)';
    ylabel 'S1_am(f)';
    grid on;

    S2_c=cos(2*pi*(fc2/fs2)*t2);
    S2_modulated1=modul(S2,S2_c,A2);

    subplot(4,3,7);
    plot(t2/fs2,S2_modulated1);
    title("Modulated S2(t)");
    xlabel 'Time(sec)';
    ylabel 'Amplitude';
    grid on;

    subplot(4,3,8);
    S2_modulated1_f=fftshift(fft(S2_modulated1));
    n2=length(S2_modulated1_f);
    freq1=(-fs2/2):(fs2/(n2-1)):(fs2/2);
    plot(freq1,(abs(S2_modulated1_f)/n2));
    title("Frequency response of modulated S2(t)");
    xlabel 'Frequency(Hz)';
    ylabel 'S2_am(f)';
    grid on;
    
    resultant_signal=S1_modulated1+S2_modulated1;
end
  
% --------------------------AM ---------------------%
function modulated = modul(message,carrier,A)
    message1=message.';
    modulated1=message1.*carrier;
    modulated=modulated1+A*carrier;
end

% --------------------------Demodulation---------------------%
function [demodulated,demodulated_filtered] = demod(modulated,fc,fm,fs,t,A)
    carrier= 2*cos(2*pi*(fc/fs)*t);
    demodulated=modulated.*carrier;
    demodulated_filtered=lowpass(demodulated,fm,fs);
    demodulated_filtered_1=demodulated_filtered-(A);
    gauss_filter=fspecial('gaussian',[1 198450],10);
    demodulated_filtered=imfilter(demodulated_filtered_1,gauss_filter);
end
