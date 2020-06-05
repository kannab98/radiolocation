function r_xx = corr(signal,m)
    r_xx = 0;
    N = length(signal);
    if m >= 0
        for n=1:N-m
            r_xx = r_xx + signal(n+m)*conj(signal(n));
        end
    else
        for n=1:N-abs(m)
            r_xx = r_xx + conj(signal(n+abs(m)))*(signal(n));
        end
    end
    r_xx = r_xx/(N - abs(m));
end

