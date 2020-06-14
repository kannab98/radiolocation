%% ������ � ��������� ������
A0 = 1;
alpha0 = -0;
x = 1 * crand(560,1);
t = zeros(length(x),1);
T = 0.00125;

for i=1:length(x)
    t(i,1) = (i-1)*T;
end

noise = x;
f0 = 10;
count_waves = 4

for j=1:length(x)
    for n=1:count_waves
        x(j) = x(j) + A0*exp(1i*2*pi*f0*n*T*(j-1)) * exp(alpha0*(j-1)*T); 
    end
end

N = length(x);
p = 64
%% ���������� ���������������� ������� ����������
p = 2*p;
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

row = zeros(p+1,1);
col = zeros(p,1);
for i=1:p
    col(i) = x(p+i);
end
for i=1:p+1
    row(i) = x(p-i+2);
end
R = [TT;TTJ]'*[TT;TTJ];

[U,S,V] = svd(R);
% % m = p;
m = count_waves;
S = diag(S);
R = zeros(p+1,p+1);
for i=1:m
    R = R + S(i)*U(:,i)*V(:,i)';
end

g = linsolve([R(:,1:p),R(:,(p+2):end)], -R(:,p+1));
a = [g(1:p);1;g(p+1:end)];
z = roots(a);

%% ����� ����������� �������� h
Z = zeros(N,length(z));
for i=1:N
    for j=1:length(z)
            Z(i,j) = z(j)^(i-1);
    end
end
h = linsolve(Z'*Z,Z'*x);

%% ���������� ���������� ���������
A = abs(h);
phi = atan(imag(h)./real(h));
f = 1/(2*pi*T) .* atan(imag(z)./real(z));
alpha = 1/T*log(abs(z));


%% �������������� �������
x1 = zeros(length(x),1);
for j=1:length(x)
    for n=1:length(z)
        x1(j) = x1(j) + A(n)*exp(1i*2*pi*f(n)*T*(j-1) + 1i*phi(n)) * exp(alpha(n)*(j-1)*T);
    end
end

x2 = zeros(length(x),1);
[A1,arg] = sort(A);
A1 = flip(A1);
arg=flip(arg);
f1 = f(arg);

phi1 = phi(arg);
alpha1 = alpha(arg);
for j=1:length(x)
    for n=1:m
        x2(j) = x2(j) + A1(n)*exp(1i*2*pi*f1(n)*T*(j-1) + 1i*phi1(n)) * exp(alpha1(n)*(j-1)*T);
    end
end 
%% �������������� �������
F = -1/(2*T):0.01:1/(2*T);
Spectrum = spectrum(F,T,p,alpha,h,f);

%% �������
figure(1); clf(); hold on;
    plot(real(x),'-')
    plot(real(x1))
    legend('������', '�����')
    hold off
 
figure(2); clf(); hold on;
    plot(real(x-noise),'-')
    plot(real(x2))
    legend('������', '�����')
    hold off
    
figure(3); clf(); hold on;
    plot(F,Spectrum)
    xlabel('f')
    ylabel('S')
%     ylim([-70,1])
hold off

figure(4); clf(); hold on;
    title('����������� ���������')
    scatter(real(z), imag(z))
    xx = -1:0.01:1;
    rho = sqrt(1-xx.^2);
    plot(xx,rho,'r')
    plot(xx,-rho,'r')
    xlabel('Re(z)')
    ylabel('Im(z)')
hold off;