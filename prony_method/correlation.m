function R = correlation(signal,p)
    R = zeros(p,p);
    N=length(signal);
    x = signal;
    for i=1:p
        for j=1:p
            for n=1:N-p
                if i<j
                    R(i,j) = R(i,j) + 1/(N-p) * ( conj(x(n+p)) * x(n) );
                else 
                    R(i,j) = R(i,j) + 1/(N-p) * ( x(n+p) * conj(x(n)) );
            
                end
            end
        end
    end
end

