x = test_function();
T = 1;

A = 0.4;
% x = A * rand(128,1);
% 
% for j=1:length(x)
%     x(j) = x(j) + exp(-2i*(j-1)*T/4) + exp(-2i*(j-1)*T/8) + exp(1i*(j-1)*T/16); 
% end

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
for i=1:length(x)/2
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

figure(1); clf(); hold on;
    plot(real(x))
    plot(real(x_model));
    legend('real','model')
hold off;


