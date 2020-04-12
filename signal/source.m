function U = source(U0,SNR,omega,t)
    U = U0/SNR*rand(1,length(t));
    for i=1:length(t)
        U(i) = U(i) + U0 * exp(1i*omega*t(i));
    end
end
