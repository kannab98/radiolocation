% Сигнал с симметричной нелинейной частотной модуляцией
clc, clear;
c=299792458;
Mn=50;
N_F1 = round(Mn*242);       % Количество отсчётов ДЗ
N_FFT_F1 = round(Mn*1*2048);
F_F1 =  3.0e6;             % Полоса частот сигнала ДЗ, Гц
Fd = Mn*473e6/156*2;         % Частота дискретизации, Гц
T_F1 =  N_F1/Fd;            % Длительность сигнала ДЗ, с
dT = 1/Fd;                  % Длительность кванта, с
%=== моделирование изменения частоты сигнала ДЗ
t1 = -(N_F1/2-1)*dT : dT : (N_F1/2-1)*dT;
Tmax = (N_F1/2-1)*dT;
y = 120000;
y1 = 220000;
y2 = 300000;
y3 = 401000;
a = .6;
b = 0.5/4*(1-a)*F_F1/((sinh(y*Tmax)));
b1 = 0.5/4*(1-a)*F_F1/((sinh(y1*Tmax)));
b2 = 0.5/4*(1-a)*F_F1/((sinh(y2*Tmax)));
b3 = 0.5/4*(1-a)*F_F1/((sinh(y2*Tmax)));
delta = a*F_F1/T_F1.*t1;
zigzag = b*(sinh(y.*t1));
zigzag1 = b1*(sinh(y1.*t1));
zigzag2 = b2*(sinh(y2.*t1));
zigzag3 = b3*(sinh(y3.*t1));
freq_F1 = delta + zigzag + zigzag1 + zigzag2 + zigzag3;
figure(1), plot(t1, freq_F1), grid on;
phase_F1 = zeros(1,N_F1);
phase_Filt = zeros(1,N_F1);
Fdop = 0000;
for k = 2:1:N_F1
    phase_Filt(1,k) = phase_Filt(1,k-1) + 2*pi*(freq_F1(1,k-1))*dT;
    phase_F1(1,k) = phase_F1(1,k-1) + 2*pi*(freq_F1(1,k-1)+Fdop)*dT;
end
win = ones(1, N_F1);
sample_F1(1,:) = win.*exp(1i*phase_F1);
sample_Filt = exp(1i*phase_Filt);
figure(2), plot(1:N_F1, real(sample_F1), 1:N_F1, imag(sample_F1)), grid on;
sample_shift_F1 = complex(zeros(1, N_FFT_F1));
sample_shift_F1(floor(N_FFT_F1/2)+1 : floor(N_FFT_F1/2)+N_F1) = sample_Filt;
S_sample_F1 = conj(fft(sample_shift_F1, N_FFT_F1));
Sx = fft(sample_F1,N_FFT_F1);
Px = Sx.*S_sample_F1;
s_compr(1,:) = (ifft(Px,N_FFT_F1));
save signal sample_F1 S_sample_F1 Mn N_F1 Fd
t2 = (-(N_FFT_F1/2)*dT : dT : (N_FFT_F1/2-1)*dT)*c/2;
figure(4),plot(t2/1e3, 20*log10((abs(s_compr))/(N_F1))); grid on;
ylim([-42 0]);
xlabel('Длительность, км');
ylabel('A, дБ');
figure(5),plot(t2/1e3, 20*log10((abs(s_compr))/(N_F1))); grid on;
ylim([-55 0]);
xlabel('Длительность, км');
ylabel('A, дБ');