addpath('./scripts');
tau = 100e-6; % 100 мкс  
df = 3e6;     % 3 МГц
step = (2*df)^-1;
t = -tau/2:step:tau/2;

figure(1); clf(); 
hold on;
    a = 0;
    plot(t,omega(t,df,tau,a));
    a = 0.5;
    plot(t,omega(t,df,tau,a));
hold off;

signal = exp(1i*omega(t,df,tau,a).* t );
figure(2); clf(); plot(t,signal);

spectrum = fft(signal);
conv_time = conv(signal, conj(signal));
conv_freq = ifft(spectrum .* conj(spectrum));

conv_time = 20*log10(abs(conv_time));
conv_freq = 20*log10(abs(conv_freq));
figure(3); clf(); 
hold on;
     plot(conv_time);
     plot(conv_freq);
     legend('t','f')
hold off


