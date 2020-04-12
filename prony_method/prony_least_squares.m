% x = test_function();
m = 1;
T = 1/100; 
A = 0.;
B = 0.01;
C = 4;
x = A * rand(64,1);
noise = x;

for j=1:length(x)
    for n=1:m
        x(j) = x(j) + exp((-B+1i)*2*pi*(j-1)*C*T*n);
    end
end

N = length(x);

%%%%%%% Корреляция
m=0;
p=4;
R = correlation(x,p);
a = auto_regressive(R)


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = length(x);
% p = N/4;
% A = ones(p);
% B = ones(p,1);
% 
% for i = 1:p
%     B(i) = 1;
%     for j = 1:p
%         A(i,j) = x(p + i - j);
%     end
% end
% 
% a = linsolve(A,B);
% a = [1; a];
% z = roots(a);
% 
% 
% alpha = log(abs(z))/T;
% f = atan( imag(z) ./ real(z) ) / (2 * pi * T);
%  
%  
% Z = zeros(p-1, p);
% for i=1:N
%     for j=1:p
%         Z(i,j) = z(j)^(i-1); 
%     end
% end
% X = x;
% h = linsolve(conj(Z')*Z,conj(Z')*X);
% 
% A = abs(h);
% theta = atan(imag(h) ./ real(h) );
% final = length(x);
% 
% x_model = zeros(length(x),1);
% for i=1:length(x)
%     for j=1:p
%         x_model(i) = x_model(i) +                 ...
%         A(j)*exp                                  ...
%             (                                     ...
%                 alpha(j)*     (i-1)*T +           ...
%                 1i*2*pi*f(j)* (i-1)*T +           ... 
%                 1i*theta(j)                       ...
%             );
%     end
% end
% 
% F = -1/(2*T):0.01:1/(2*T);
% S = spectrum(F,T,p,alpha,h,f);
% figure(1); clf(); hold on;
%     plot(real(x(1:end/2)));
%     plot(real(x_model(1:end/2)));
%     plot(real(noise(1:end/2)));
%     legend('real signal','Prony`s method','noise')
% hold off;
% 
% figure(2); clf(); hold on;
%     plot(F,S)
%     xlim([-1,+C*N+2])
% hold off;
