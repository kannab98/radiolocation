% Some data - replace it with yours (its from an earlier project)
x = -1:0.01:5;
y = exp(x); 




f = 3e6;
t = 0:1e-8:2e-6;
A = 0.01;
x = A * rand(length(t),1);
for j=1:length(x)
        x(j) = x(j) + exp(1j*2*pi*f*t(j));
end
x=x;
t=t';
plot(t,x)

% 
% x=x';
N =1;
% 
% args=[];
Nargs = 2;
% 
arg=cell(Nargs,N);
alphabet=char(97:122);
command = [];
for i=1:N
    for j=1:Nargs
        args{i,j}=[alphabet(j),num2str(i)];
    end
end

for i=1:N
    command=[command,[args{i,1},'*exp(',args{i,2},'*x*1j)+']];
end

command = command(1:end-1);
fitfun = fittype(command)
p0 = [1,1];
low = [0,0];
up = [1,1];
[fitted_curve,gof] = fit(t,x,fitfun,'StartPoint',p0,'Lower',low)

% Save the coeffiecient values for a,b,c and d in a vector
coeffvals = coeffvalues(fitted_curve);
% Plot results
scatter(t, x, 'r.')
hold on
plot(x,fitted_curve(x))
hold off