function signal = receiver(N,U,lambda,d,h,alpha,earth)
    signal = zeros(N,length(U));
    for i=1:N
        for j=1:length(U)        
            signal(i,j) = U(j) .* exp(1i*2*pi/lambda * d * sin(alpha) * (i-1) ); 
        end
    end

    if earth == 1
        signal1 = zeros(N,length(U));
        for i=1:N
            for j=1:length(U)        
                signal1(i,j) = U(j) .* exp( 1i*2*pi/lambda * ....
                                            (2*h + d*(i-1) ) * sin(alpha)...
                                         ); 
            end
        end
        signal = cat(1,signal,signal1);
    end
end

