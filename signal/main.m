set(0,'DefaultTextInterpreter',          'latex');
set(0,'DefaultLegendInterpreter',        'latex');
set(0,'DefaultAxesTickLabelInterpreter', 'latex');
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultTextFontSize',12);

%% Параметры РЛС

N = 4; % Число строк в решетке
d = 1e-1; % Расстояние между строками в метрах.
h = 0.5; % Высота приемника над поверхностью земли в метрах.

%% Параметры источника

U0 = 1; % Амплитуда сигнала в у.е. (В, дБ?).
f = 3e8; % Частота в Гц.
omega = 2*pi*f; % Частота сигнала в 1/рад.
c = 3e8; % Скорость света в вакууме в м/с.
lambda = omega/(2*pi*c); % Длина волны сигнала в метрах, среда без дисперсии.
t = 0:0.001:4; % Время измерения в долях от периода сигнала
T = 1/f;
t = t .* T; 
alpha = pi/6; % Теоретический угол места источника. 
SNR = 10e6;
U = source(U0,SNR,omega,t);
%% Моделирование принятого сигнала на РЛС;
% Переменная earth включает или отключает отраженный сигнал.
earth = 0;
signal = receiver(N,U,lambda,d,h,alpha,earth); % Сигнал от источника
figure(1); clf(); hold on;
for i=1:N
    plot(t,real(signal(i,:)),'DisplayName',['$R^0_',num2str(i),'$'] )
    xlabel('$t, s$')
    ylabel('$U, V$')
end
if earth == 1
    for i=N+1:2*N
        plot(t+ 2*h*sin(alpha)/lambda/omega , real(signal(i,:)),'--','DisplayName',['$R^1_',num2str(i),'$'])
        xlabel('$t, s$')
        ylabel('$U, V$')
    end
end
legend();
hold off;


%% Вычисление угла места;
% В функции elevation_calc по идее также должна выполняться спектральная
% оценка сигнала signal. То есть аргумент lambda должен искаться внутри самой 
% функции.
% Пока что предполагаю спектр известным.

elevation = elevation_calc(signal,lambda,d)
