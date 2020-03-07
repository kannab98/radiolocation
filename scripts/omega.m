function omega = omega(t,df,tau,a)  
    omega =  2*pi*df * (t/tau + a*sinh(t/tau));
end