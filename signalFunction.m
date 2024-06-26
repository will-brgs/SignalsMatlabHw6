function [r, SNR_arr, xn] = signalFunction(bitrate, sigma_arr)
    
    Tp = 0.1; % Half pulse width
    sample_period = Tp/50; % dt, pulse and recieve sample period
    sample_freq = 1/sample_period; % Frequency of pulse and recieve signal 
    
    bit_rate = bitrate; %Fb, frequency of bits sent out
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
    r = zeros(length(y), length(sigma_arr));
    % figure, subplot (2,1,1),plot(y)
    % subplot(2,1,2),stem(xn)
    SNR_arr = zeros(1, length(sigma_arr));
    for i = 1:length(sigma_arr)
        sigma = sigma_arr(i);
        noise = sigma * max(y) * randn(1,length(y));
        r(:, i) = y + (noise * sigma);
        SNR_arr(i) = (sum(y.^2))/(sum(noise.^2));
    end
    
    SNR_arr(1) = 1;

end