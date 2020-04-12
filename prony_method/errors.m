function e = errors(signal,p)
R = correlation(signal,p);
N = length(signal);
rho=zeros(N,1);
k=zeros(N,1);
a=zeros(N,1);
rho(1)= R(1,1);
for i=2:p
    for j=1:p
        k(i) = k(i) + 1/rho(i-1) * (-1)  * ( a(j)*R(1,p-j) );
    end
    rho(i) = rho(i-1) * (1 - abs(k)^2);
%     a(i) = ...

errf = signal;
errb = signal;
for i=1:N
    errb(i) = signal(i-p) + k
    for j=1:p
    

for i=1:N
    errf(i) = signal(i) 
end

end




