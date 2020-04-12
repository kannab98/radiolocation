set(0,'DefaultTextInterpreter',          'latex');
set(0,'DefaultLegendInterpreter',        'latex');
set(0,'DefaultAxesTickLabelInterpreter', 'latex');
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultTextFontSize',12);

%% ��������� ���

N = 4; % ����� ����� � �������
d = 1e-1; % ���������� ����� �������� � ������.
h = 0.5; % ������ ��������� ��� ������������ ����� � ������.

%% ��������� ���������

U0 = 1; % ��������� ������� � �.�. (�, ��?).
f = 3e8; % ������� � ��.
omega = 2*pi*f; % ������� ������� � 1/���.
c = 3e8; % �������� ����� � ������� � �/�.
lambda = omega/(2*pi*c); % ����� ����� ������� � ������, ����� ��� ���������.
t = 0:0.001:4; % ����� ��������� � ����� �� ������� �������
T = 1/f;
t = t .* T; 
alpha = pi/6; % ������������� ���� ����� ���������. 
SNR = 10e6;
U = source(U0,SNR,omega,t);
%% ������������� ��������� ������� �� ���;
% ���������� earth �������� ��� ��������� ���������� ������.
earth = 0;
signal = receiver(N,U,lambda,d,h,alpha,earth); % ������ �� ���������
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


%% ���������� ���� �����;
% � ������� elevation_calc �� ���� ����� ������ ����������� ������������
% ������ ������� signal. �� ���� �������� lambda ������ �������� ������ ����� 
% �������.
% ���� ��� ����������� ������ ���������.

elevation = elevation_calc(signal,lambda,d)
