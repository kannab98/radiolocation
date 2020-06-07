 %% Сигнал и параметры модели
A0 = 1;
alpha0 = -0.2;
x = 0. * rand(560,1);
t = zeros(length(x),1);
T = 0.00125;

for i=1:length(x)
    t(i,1) = (i-1)*T;
end

noise = x;
f0 = 2;
count_waves = 2;

for j=1:length(x)
    for n=1:count_waves
        x(j) = x(j) + A0*exp(1i*2*pi*f0*n*T*(j-1)) * exp(alpha0*(j-1)*T); 
    end
end

% T = 1;
% x = test_function();
N = length(x);
p = 2;
%% Попытка учитывания известных компонент
c = fliplr(poly(z0));
y = zeros(p,1);

for n=p+1:2*p
    for k=1:q+1
        y(n-p) = y(n-p) + c(k) .* x(n-k+1);
    end
end

% % Составление тёплицевой матрицы сигнала
% row = zeros(p+1,1);
% col = zeros(p,1);
% for i=1:p
%     col(i) = x(p+i);
% end
% for i=1:p+1
%     row(i) = x(p-i+2);
% end
%     
% TT = toeplitz(col, row);
% TTJ = conj(fliplr(TT));
% 
% % Модифицированная матрица  ковариации 
% % учитывается ошибка предсказания назад и вперед
% R = [TT;TTJ]'*[TT;TTJ];
% g = linsolve([R(:,1:p),R(:,(p+2):end)], -R(:,p+1));
% a = [g(1:p);1;g(p+1:end)];
% %% Вычисление матрицы ковариации (см. формулу (6) в документе)

R_xx = corr_mat(x,p+1);
R0 = R_xx(:,1);
R_xx = R_xx(:,2:end);

%% Коэффициенты полинома (см. формулу (6) в документе)
a = linsolve(R_xx,-R0);
a = [1;a];
z = roots(a);
%% Поиск комплексных амплитуд h
Z = zeros(N,p);
for i=1:N
    for j=1:p
            Z(i,j) = z(j)^(i-1);
    end
end
h = linsolve(Z'*Z,Z'*x);


%% Вычисление параметров экспонент
A = abs(h);
phi = atan(imag(h)./real(h));
% Определяется верно
f = 1/(2*pi*T) .* atan(imag(z)./real(z));
alpha = 1/T*log(abs(z));


x1 = zeros(length(x),1);
for j=1:length(x)
    for n=1:p
        x1(j) = x1(j) + A(n)*exp(1i*2*pi*f(n)*T*(j-1) + 1i*phi(n)) * exp(alpha(n)*(j-1)*T);
    end
end

F = -1/(2*T):0.01:1/(2*T);
S1 = spectrum(F,T,p,alpha,h,f);
S2 = spectrum2(F,T,p,alpha,h,f);
figure(1); clf(); hold on;
    plot(real(x))
    plot(real(x1))
    legend('начальные данные', 'прони')
    hold off
    
figure(2); clf(); hold on;
    plot(F,S1)
    plot(F,S2)
    ylim([-10,1])
hold off

figure(3); clf(); hold on;
    scatter(real(z), imag(z))
    xx = -1:0.01:1;
    rho = sqrt(1-xx.^2);
    plot(xx,rho,'r')
    plot(xx,-rho,'r')
hold off;