function omega = omega(t,df,tau,a)
    x = t/tau;
    omega =  2*pi*df * ( x + sinh(a*x) );
end