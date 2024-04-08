%% Signals and Systems Matlab Homework #6
%% Introduction
% * Author:                   Will Burgess, 
% * Class:                    ESE 351
% * Date:                     Created 3/25/2024, Last Edited 2/2/2024
%% Housekeeping
close all
clear
clc
%% Part  1: Create Binary PAM System
%% Generate Input Signal and Add Noise Factor, Bitrate = 1/Tp
Tp = 0.1; % Half pulse width
sample_period = Tp/50; % dt, pulse and recieve sample period
sample_freq = 1/sample_period; % Frequency of pulse and recieve signal 

bit_rate = 1/(1 * Tp); %Fb, frequency of bits sent out
bit_period = 1/bit_rate; % Ts, Time between bits sent out

rect = ones(1,50);
pulse = 2 * conv(rect, rect);
pulse_fft = fftshift(pulse);

%figure, hold on
% subplot(3,1,1),stem(x);
% subplot(3,1,2),plot(n);
% subplot(3,1,3),plot(r); 
%stem(pulse);
%plot(mag2db(abs(pulse_fft)));
%hold off

N = 20;

maxTime = N * bit_period;

xn = 2 * ((rand(1, N) > 0.5) - 0.5);
a = 0;
imp_train = zeros(1,N * bit_period * sample_freq);
for k = 1:length(imp_train)
    if mod(k - 1, sample_freq * bit_period) == 0
    a = a + 1;
    imp_train(k) = xn(a);
    else
    imp_train(k) = 0;    
    end
end

sampleTimes = 0:sample_period:(N*bit_period)-sample_period;

y = conv(imp_train,pulse);
% figure, subplot (2,1,1),plot(y)
% subplot(2,1,2),stem(xn)

sigma = 1;
noise = sigma * max(y) * randn(1,length(y));
r = y + (noise * sigma);

%figure,plot(r)
%% Part f:  Create various outputs

% i, p(t), p(w)
figure, hold on
subplot(2,1,1);
stem(pulse)
ylabel('Amplitude')
xlabel('Index')
title('Pulse P(t), Bitrate = 1/Tp')



subplot(2,1,2)
plot(pulse_fft)
ylabel('Amplitude')
xlabel('Frequency (Hz)')
title('Pulse Sprectum P(w), Bitrate = 1/Tp')
hold off

times = linspace(0, maxTime + 2 * Tp, length(y));

% ii, y(t)
figure, hold on
plot(times, y)
xlim([0 maxTime + 2 * Tp]);
ylabel('Amplitude')
xlabel('Time (s)')
title('Transmitted Signal y(t), Bitrate = 1/Tp')
hold off

% iii, r(t)
figure, hold on
plot(times,r)
xlim([0 maxTime + 2 * Tp]);
ylabel('Amplitude')
xlabel('Time (s)')
title('Recieved Signal r(t), Bitrate = 1/Tp')
hold off

% iv, Sent Messave vs. decoded message
times_sent = 0:bit_period:maxTime-bit_period;
figure, hold on
subplot(2,1,1)
stem(times_sent,xn)
ylabel('Amplitude')
xlabel('Time (s)')
title('Transmitted Signal xn, Bitrate = 1/Tp')
subplot(2,1,2)

filtered = conv(r,pulse);
decoded = zeros(1, N);

a = 0;
pulselen = length(pulse);
filterlen = length(filtered);

factor = 1/(bit_rate * Tp); % find factor relating Ts and Tp, use that to modify pulselen
% please name this something other than factor

% don't ask me to explain this
for i = pulselen + 1:(pulselen * factor + mod(factor, 2))/2:filterlen-pulselen * factor - 1
    a = a + 1;
    if(filtered(i) > 0)
        decoded(a) = 1;
    else
       decoded(a) = -1;
    end
end

stem(times_sent,decoded);
ylabel('Amplitude')
xlabel('Time (s)')
title('Decoded Signal r(t), Bitrate = 1/Tp')
hold off

error = (sum(xn ~= decoded)/length(decoded)); 
SNR = (sum(y.^2))/(sum(noise.^2));

disp(['Bitrate: ' ,num2str(bit_rate), ' bits/second'])
disp(['Standard Deviation: ' , num2str(sigma)])
disp(['SNR: ' , num2str(SNR)])
disp(['Error: ' ,num2str(error),' percent'])


%% Generate Input Signal and Add Noise Factor, Bitrate = 1/2*Tp
Tp = 0.1; % Half pulse width
sample_period = Tp/50; % dt, pulse and recieve sample period
sample_freq = 1/sample_period; % Frequency of pulse and recieve signal 

bit_rate = 1/(2 * Tp); %Fb, frequency of bits sent out
bit_period = 1/bit_rate; % Ts, Time between bits sent out

rect = ones(1,50);
pulse = 2 * conv(rect, rect);
pulse_fft = fftshift(pulse);

%figure, hold on
% subplot(3,1,1),stem(x);
% subplot(3,1,2),plot(n);
% subplot(3,1,3),plot(r); 
%stem(pulse);
%plot(mag2db(abs(pulse_fft)));
%hold off

N = 20;

maxTime = N * bit_period;

xn = 2 * ((rand(1, N) > 0.5) - 0.5);
a = 0;
imp_train = zeros(1,N * bit_period * sample_freq);
for k = 1:length(imp_train)
    if mod(k - 1, sample_freq * bit_period) == 0
    a = a + 1;
    imp_train(k) = xn(a);
    else
    imp_train(k) = 0;    
    end
end

sampleTimes = 0:sample_period:(N*bit_period)-sample_period;

y = conv(imp_train,pulse);
% figure, subplot (2,1,1),plot(y)
% subplot(2,1,2),stem(xn)

sigma = 1;
noise = sigma * max(y) * randn(1,length(y));
r = y + (noise * sigma);

%figure,plot(r)
%% Part f:  Create various outputs

% i, p(t), p(w)
figure, hold on
subplot(2,1,1);
stem(pulse)
ylabel('Amplitude')
xlabel('Index')
title('Pulse P(t), Bitrate = 1/(2*Tp)')

f = (0:1:length(pulse_fft)-1);


subplot(2,1,2)
plot(f,pulse_fft)
ylabel('Amplitude')
xlabel('Frequency (Hz)')
title('Pulse Sprectum P(w), Bitrate = 1/(2*Tp)')
hold off

times = linspace(0, maxTime + 2 * Tp, length(y));

% ii, y(t)
figure, hold on
plot(times, y)
xlim([0 maxTime + 2 * Tp]);
ylabel('Amplitude')
xlabel('Time (s)')
title('Transmitted Signal y(t), Bitrate = 1/(2*Tp)')
hold off

% iii, r(t)
figure, hold on
plot(times,r)
xlim([0 maxTime + 2 * Tp]);
ylabel('Amplitude')
xlabel('Time (s)')
title('Recieved Signal r(t), Bitrate = 1/(2*Tp)')
hold off



% iv, Sent Messave vs. decoded message
times_sent = 0:bit_period:maxTime-bit_period;
figure, hold on
subplot(2,1,1)
stem(times_sent,xn)
ylabel('Amplitude')
xlabel('Time (s)')
title('Transmitted Signal xn, Bitrate = 1/(2*Tp)')
subplot(2,1,2)

filtered = conv(r,pulse);
decoded = zeros(1, N);

a = 0;
pulselen = length(pulse);
filterlen = length(filtered);

factor = 1/(bit_rate * Tp); % find factor relating Ts and Tp, use that to modify pulselen
% please name this something other than factor

% don't ask me to explain this
for i = pulselen + 1:(pulselen * factor + mod(factor, 2))/2:filterlen-pulselen * factor - 1
    a = a + 1;
    if(filtered(i) > 0)
        decoded(a) = 1;
    else
       decoded(a) = -1;
    end
end

stem(times_sent,decoded);
ylabel('Amplitude')
xlabel('Time (s)')
title('Decoded Signal r(t), Bitrate = 1/(2*Tp)')
hold off

error = (sum(xn ~= decoded)/length(decoded)); 
SNR = (sum(y.^2))/(sum(noise.^2));

disp(['Bitrate: ' ,num2str(bit_rate*100), ' bits/second'])
disp(['Standard Deviation: ' , num2str(sigma)])
disp(['SNR: ' , num2str(SNR)])
disp(['Error: ' ,num2str(error),' percent'])

%% Part 2: Performance Test
sigma_arr = [0, 1, 2, 3, 4];
r_Tp = signalFunction(1/Tp, sigma_arr);
r_2Tp = signalFunction(1/(2*Tp), sigma_arr);
times_sent_Tp = 0:Tp:N*Tp-Tp;
figure, hold on
for j = 1:5
    filtered_Tp = conv(r_Tp(:, j),pulse);
    a = 0;
    decoded_Tp = zeros(1, N);
    pulselen = length(pulse);
    filterlen = length(filtered_Tp);
    
    factor = 1/((1/Tp) * Tp); % find factor relating Ts and Tp, use that to modify pulselen
    % please name this something other than factor
    
    % don't ask me to explain this
    for i = pulselen + 1:(pulselen * factor + mod(factor, 2))/2:filterlen-pulselen * factor - 1
        a = a + 1;
        if(filtered_Tp(i) > 0)
            decoded_Tp(a) = 1;
        else
           decoded_Tp(a) = -1;
        end
    end
    
    
    subplot(3, 2, j);
    stem(times_sent_Tp, decoded_Tp);
    title(['Sigma = ',num2str(sigma_arr(j))])
    sgtitle('Matched Fitler With Bitrate=1/Tp')
    xlabel('Time (s)')
    ylabel('Decoded Output')
end
hold off;


times_sent_2Tp = 0:2*Tp:2*N*Tp-(2*Tp);
figure, hold on
for j = 1:5
    filtered_2Tp = conv(r_2Tp(:, j),pulse);
    a = 0;
    decoded_2Tp = zeros(1, N);
    pulselen = length(pulse);
    filterlen = length(filtered_2Tp);
    
    factor = 1/((1/(2*Tp)) * Tp); % find factor relating Ts and Tp, use that to modify pulselen
    % please name this something other than factor
    
    % don't ask me to explain this
    for i = pulselen + 1:(pulselen * factor + mod(factor, 2))/2:filterlen-pulselen * factor - 1
        a = a + 1;
        if(filtered_2Tp(i) > 0)
            decoded_2Tp(a) = 1;
        else
           decoded_2Tp(a) = -1;
        end
    end
    
    
    subplot(3, 2, j);
    stem(times_sent_2Tp, decoded_2Tp);
    title(['Sigma = ',num2str(sigma_arr(j))])
    sgtitle('Matched Fitler With Bitrate=1/(2*Tp)')
    xlabel('Time (s)')
    ylabel('Decoded Output')
end
hold off

figure, hold on
for j = 1:5
    a = 0;
    unfiltered_2Tp = zeros(1, N);
    pulselen = length(pulse);
    filterlen = length(r_2Tp(:,j));
    
    factor = 1/((1/(2*Tp)) * Tp); % find factor relating Ts and Tp, use that to modify pulselen
    % please name this something other than factor
    
    % don't ask me to explain this
    for i = pulselen + 1:(pulselen * factor + mod(factor, 2))/2:filterlen-pulselen * factor - 1
        a = a + 1;
        if(r_2Tp(i) > 0)
            unfiltered_2Tp(a) = 1;
        else
           unfiltered_2Tp(a) = -1;
        end
    end
    
    
    subplot(3, 2, j);
    stem(times_sent_2Tp, unfiltered_2Tp);
    title(['Sigma = ',num2str(sigma_arr(j))])
    sgtitle('Sign-based Reciever With Bitrate=1/(2*Tp)')
    xlabel('Time (s)')
    ylabel('Decoded Output')
end
hold off
