tau = 100e-6; % 100 мкс  
df = 3e6;     % 3 МГц
step = (2*df)^-1;
t = 0:step:tau;

y = sin(omega(t,df,tau) .* t );
plot(t,y)

spectrum = fft(y);
Y1 = conv(y,y);
Y1 = Y1(1:length(y));
Y2 = ifft(spectrum .* spectrum);

figure(2);
clf();
title('Функция корреляции')
t1=tau:-step:0;
hold on;
plot(t1,Y1);
plot(t1,Y2);
legend('t-способ','f-способ')
hold off;


function omega = omega(t,df,tau)
    omega =  2*pi*df/tau * t;
end