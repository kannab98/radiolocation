 %% Сигнал и параметры модели
A0 = 1;
alpha0 = -0;
x = 0.0 * rand(256,1);
t = zeros(length(x),1);
T = 0.00125;

for i=1:length(x)
    t(i,1) = (i-1)*T;
end

noise = x;
f0 = 2;
count_waves = 4;

for j=1:length(x)
    for n=1:count_waves
        x(j) = x(j) + A0*exp(1i*2*pi*f0*n*T*(j-1)) * exp(alpha0*(j-1)*T); 
    end
end
%%
% T = 1;
% x = test_function();
N = length(x);
p = N/2;


% Составление тёплицевой матрицы сигнала
row = zeros(p+1,1);
col = zeros(p,1);
for i=1:p
    col(i) = x(p+i);
end
for i=1:p+1
    row(i) = x(p-i+2);
end
    
TT = toeplitz(col, row);
TTJ = conj(fliplr(TT));

% Модифицированная матрица  ковариации 
% учитывается ошибка предсказания назад и вперед
% R = [TT;TTJ]'*[TT;TTJ];
% g = linsolve([R(:,1:p),R(:,(p+2):end)], -R(:,p+1));
% a = [g(1:p);1;g(p+1:end)];

% Обычный метод, только ошибка вперед
a1 = linsolve(TT(:,2:end), - TT(:,1));
a = [1; a1];



% Р­С‚Р°Рї 2 
z = roots(a);

alpha = log(abs(z))/T;
f = atan( imag(z) ./ real(z) ) / (2 * pi * T);
 

% Р­С‚Р°Рї 3  
% p = 2*p-1;
Z = zeros(p-1, p);
for i=1:p
    for j=1:p
        Z(i,j) = z(j)^(i-1); 
    end
end
X = x(1:p);
h = linsolve(Z,X);

A = abs(h);
theta = atan(imag(h) ./ real(h) );
final = length(x);

x_model = zeros(length(x),1);
for i=1:length(x)
    for j=1:p
        x_model(i) = x_model(i) +                 ...
        A(j)*exp                                  ...
            (                                     ...
                alpha(j)*     (i-1)*T +           ...
                1i*2*pi*f(j)* (i-1)*T +           ... 
                1i*theta(j)                       ...
            );
    end
end

F = -1/(2*T):0.01:1/(2*T);
S1 = spectrum(F,T,p,alpha,h,f);
S2 = spectrum2(F,T,p,alpha,h,f);
figure(1); clf(); hold on;
    plot(real(x(1:end/2)));
    plot(real(x_model(1:end/2)));
    plot(real(noise(1:end/2)));
    legend('real signal','Prony`s method','noise')
hold off;

figure(2); clf(); hold on;
    plot(F,S1)
    plot(F,S2)
    ylim([-30,1])
hold off;
