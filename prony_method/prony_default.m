  % x = test_function();
N = 4;
T = 1/800; 
A = 0.0;
B = 0.0;
C = 4;
x = A * rand(256,1);
noise = x;
for j=1:length(x)
    for n=1:N
        x(j) = x(j) + exp((-B/n+1i)*2*pi*(j-1)*C*T*n);
    end
end

p = length(x)/2;
A = ones(p);
B = ones(p,1);

for i = 1:p
    B(i) = - x(p+i);
    for j = 1:p
        A(i,j) = x(p + i - j);
    end
end

a = linsolve(A,B);
a = [1; a];
z = roots(a);

alpha = log(abs(z))/T;
f = atan( imag(z) ./ real(z) ) / (2 * pi * T);
 
 
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
S = spectrum(F,T,p,alpha,h,f);
figure(1); clf(); hold on;
    plot(real(x(1:end/2)));
    plot(real(x_model(1:end/2)));
    plot(real(noise(1:end/2)));
    legend('real signal','Prony`s method','noise')
hold off;

figure(2); clf(); hold on;
    plot(F,S)
    xlim([0,+C*N+2])
%   ylim([0,max(S)+1])
hold off;