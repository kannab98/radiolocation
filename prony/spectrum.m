function S = spectrum(F,T,p,alpha,h,f)
    N = length(F);
    X = zeros(1,N);

    z = @(f,alpha) exp( alpha.*T + 1i.*2.*pi.*f.*T);
    
    for i=1:N
        for j=1:length(f)
            z0 = z(F(i),0);
            zk = z(f(j),alpha(j));
            if abs(zk/z0 - 1) < 1e-6 
                X(i) = X(i) + 1e6;
            else
                X(i) = X(i) + h(j)/(1 - zk/z0 ) - h(j)/(1 - 1/conj(zk) * 1/z0 );
            end
        end
    end
    S = abs(X) .^ 2;
    S = 10 .* log10(S/max(S));
end

