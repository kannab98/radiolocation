function r_xx = corr(signal,m)
    r_xx = 0;
    N = length(signal);
    if m >= 0
        for n=1:N-m-1
            r_xx = r_xx + signal(n+m+1)*conj(signal(n+1));
        end
    else
        for n=1:N-abs(m)-1
            r_xx = r_xx + conj(signal(n+abs(m)+1))*(signal(n+1));
        end
    end
    r_xx = r_xx/(length(signal) - abs(m) -1);
end

