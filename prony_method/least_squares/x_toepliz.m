function  X = x_toepliz(signal, p)
    N = length(signal);
    X = zeros(N+p,p+1);
    
    for i=1:N+p
        for j=1:p+1
            if (1+i-j < 1) || (1+i-j > N)
                X(i,j) = 0;
            else
                X(i,j) = signal(1+i-j);  
            end
        end
    end
    
end

