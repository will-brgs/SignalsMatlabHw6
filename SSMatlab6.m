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
%% Generate Input Signal and Add Noise Factor
Tp = 0.1; % Half pulse width
sample_period = Tp/50; % dt, pulse and recieve sample period
sample_freq = 1/sample_period; % Frequency of pulse and recieve signal 

bit_rate = 1/(Tp); %Fb, frequency of bits sent out
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
title('Pulse P(t)')

subplot(2,1,2)
plot(pulse_fft)
ylabel('Amplitude')
xlabel('Frequency (rad)')
title('Pulse Sprectum P(w)')
hold off

times = linspace(0, maxTime + 2 * Tp, length(y));

% ii, y(t)
figure, hold on
plot(y)
ylabel('Amplitude')
xlabel('Time (s)')
title('Transmitted Signal y(t)')
hold off

% iii, r(t)
figure, hold on
plot(times,r)
ylabel('Amplitude')
xlabel('Time (s)')
title('Recieved Signal r(t)')
hold off

% iv, Sent Messave vs. decoded message
times_sent = 0:bit_period:2-bit_period;
figure, hold on
subplot(2,1,1)
stem(times_sent,xn)
ylabel('Amplitude')
xlabel('Time (s)')
title('Transmitted Signal xn')
subplot(2,1,2)

filtered = conv(r,pulse);

decoded = zeros(length(N));
a = 0;
pulselen= length(pulse)
filterlen= length(filtered)


for i = 100:50:length(filtered)-100
    if(filtered(i) > 0)
        a = a + 1;
    decoded(a) = 1;
    else
       a = a + 1;
       decoded(a) = -1;
    end
end

stem(times_sent,decoded);
ylabel('Amplitude')
xlabel('Time (s)')
title('Decoded Signal r(t)')
hold off
%% Part 2: Performance Test
