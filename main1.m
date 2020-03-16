x = test_function();
T = 1;
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
 
% 
Z = zeros(p-1, p);
for i=1:p
    for j=1:p
        Z(i,j) = z(i)^(j-1); 
    end
end
X = x(1:p);
h = linsolve(Z,X);

A = abs(h);
theta = atan(imag(h) ./ real(h) );
