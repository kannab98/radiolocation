addpath('./scripts');
set(0,'DefaultTextInterpreter',          'latex');
set(0,'DefaultLegendInterpreter',        'latex');
set(0,'DefaultAxesTickLabelInterpreter', 'latex');
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultTextFontSize',12);

%begin{constants}
    tau = 100e-6; % 100 microseconds
    df = 3e6;     % 3 MHz
    freq_doppler = df/100;
    step = (200*df)^-1;
    t = -tau/2:step:tau/2;
    c = 3e8; % light   velocity, m/s
    v = 10;   % doppler velocity, m/s
%end{constants}

figure(1); clf(); 
hold on; grid on; grid minor;
    a = 0;
    plot(t,omega(t,df,tau,a));
    a = 2;
    plot(t,omega(t,df,tau,a));
    xlabel('$t$')
    ylabel('$\omega(t)$')
    legend('$a$','$b$')
    savepdf('fig/omega.pdf')
hold off;

%begin{linear}
a=0;
signal = exp(1i .* omega(t,df,tau,a) .* t );
signal_doppler = exp(1i .* omega(t,df,tau,a) .* t + 1i .* freq_doppler .* t);
% phase = zeros(1,length(signal));
% for i=2:length(signal)
%    phase(i) = 1i .* ( omega(t(i),df,tau,a) + freq_doppler ) .* t(i) ;
% end
% signal_doppler = exp(phase);
figure(2); clf(); hold on;
plot(t,signal);
plot(t,signal_doppler);
hold off;

ylim([0.2,1.05]);
xlim([-tau/8,tau/8]);
xlabel('$t,s$');
ylabel('amplitude');
legend('$a$','$b$')
savepdf('fig/signal.pdf');


spectrum = fft(signal);
spectrum_doppler = fft(signal_doppler);

conv_time = conv(signal, conj(signal));
conv_freq = ifft(spectrum .* conj(spectrum));
conv_time_doppler = conv(signal, conj(signal_doppler));
conv_freq_doppler = ifft(spectrum .* conj(spectrum_doppler));

conv_time = 20*log10(abs(conv_time));
conv_freq = 20*log10(abs(conv_freq));
conv_time_doppler = 20*log10(abs(conv_time_doppler));
conv_freq_doppler = 20*log10(abs(conv_freq_doppler));

figure(3); clf(); hold on; grid on; grid minor;
    t1=-tau/2:step:tau/2 + tau;
    plot(t1,conv_time);
%     plot(t,conv_freq);
    plot(t1,conv_time_doppler);
%     plot(t,conv_freq_doppler);
    legend('$t$','$f$')
    xlabel('$t, s$')
    ylabel('$x(t)*x(t),~dB$')
    legend('$a$','$b$')
    ylim([65,max(conv_time)])
    xlim([4.96e-5,5.04e-5])
    savepdf('fig/linear_conv.pdf')
hold off

%end{linear}


%begin{nonlinear}
a=2;
signal = exp(1i .* omega(t,df,tau,a) .* t );
signal_doppler = exp(1i .* omega(t,df,tau,a) .* t + 1i .* freq_doppler .* t);
spectrum = fft(signal);
spectrum_doppler = fft(signal_doppler);

conv_time = conv(signal, conj(signal));
conv_freq = ifft(spectrum .* conj(spectrum));
conv_time_doppler = conv(signal, conj(signal_doppler));
conv_freq_doppler = ifft(spectrum .* conj(spectrum_doppler));

conv_time = 20*log10(abs(conv_time));
conv_freq = 20*log10(abs(conv_freq));
conv_time_doppler = 20*log10(abs(conv_time_doppler));
conv_freq_doppler = 20*log10(abs(conv_freq_doppler));

figure(4); clf(); hold on; grid on; grid minor;
    t1=-tau/2:step:tau/2 + tau;
    plot(t1,conv_time);
%     plot(t,conv_freq);
    plot(t1,conv_time_doppler);
%     legend('$t$','$f$')
    xlabel('$t, s$')
    ylabel('$x(t)*x(t),~dB$')
    legend('$a$','$b$')
    ylim([65,max(conv_time)])
    xlim([4.98e-5,5.02e-5])
    savepdf('fig/nonlinear_conv.pdf')
hold off
%end{nonlinerar}