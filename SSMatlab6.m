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
Tp = 0.1;
dt = Tp/50;

Fs = 1/dt;
maxTime = Tp;
sampleTimes = 0:1/Fs:maxTime;

x = round(rand(1,length(sampleTimes)));
x(x==0) = -1;
n = randn(1,length(x));

sigma = 1;

n = n * sigma;

r = x + n;

rect = ones(1,50);

pulse = 2 * conv(rect, rect);

pulse_fft = fftshift(pulse);

figure, hold on
% subplot(3,1,1),stem(x);
% subplot(3,1,2),plot(n);
% subplot(3,1,3),plot(r); 
%stem(pulse);
plot(mag2db(abs(pulse_fft)));
hold off

N = 100;
bits = 2 * ((rand(1, N) > 0.5) - 0.5);

%% Part 2: Performance Test