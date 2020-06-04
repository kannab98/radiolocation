 %% Сигнал и параметры модели
A0 = 0.0;
alpha0 = 0.001;
x = A0 * rand(160,1);
t = zeros(length(x),1);
T = 0.0125;
for i=1:length(x)
    t(i,1) = (i-1)*T;
end

noise = x;
f0 = 2;
count_waves = 1;

for j=1:length(x)
    for n=1:count_waves
        x(j) = x(j) + exp(1i*2*pi*f0*n*T*(j-1)) * exp(alpha0*(j-1)*T); 
    end
end

N = length(x);
p = count_waves ;
%% Вычисление матрицы ковариации (см. формулу (6) в документе)
R_xx = corr_mat(x,p+1);
R0 = R_xx(:,1);
R_xx = R_xx(:,2:end);

%% Коэффициенты полинома (см. формулу (6) в документе)
a = linsolve(R_xx,R0);
a = [1;a];
z = roots(a);
%% Поиск комплексных амплитуд h
Z = zeros(p,p);

% for i=1:p
%     for j=1:p
%         for n=1:N
%             Z(i,j) = Z(i,j) + z(i)^(n-1)*conj(z(j))^(n-1);
%         end
%     end
% end
% 
% C = zeros(p,1);
% 
% for i=1:p
%     for n=1:N
%         C(i) = C(i) + x(n)*conj(z(i))^(n-1);
%     end
% end
% h = linsolve(Z,C)

% Более короткая запись, по сравнению с блоком выше

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
    for n=1:count_waves
        x1(j) = x1(j) + A(n)*exp(1i*2*pi*f(n)*T*(j-1) + 1i*phi(n)) * exp(alpha0*(j-1)*T);
    end
end

scatter(t, real(x), 'r+')

plot(t, real(x))
hold on
plot(t, real(x1))
legend('начальные данные', 'прони')
hold off