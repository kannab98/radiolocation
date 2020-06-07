function S = spectrum(F,T,p,alpha,h,f)
    N = length(F);
    X = zeros(1,N);

    z = @(f,alpha) exp( alpha.*T + 1i.*2.*pi.*f.*T);
    
    for i=1:N
        for j=1:p
            X(i) = X(i) + h(j)/(1 - z(f(j),alpha(j))/z(F(i),0) );
        end
    end
    S = abs(T .* X) .^ 2;
    S = 10 .* log10(S/max(S));
end

