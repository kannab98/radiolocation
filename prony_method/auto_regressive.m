function a = auto_regressive(correlation_matrix)
    rho = zeros(length(correlation_matrix),1);
    rho(1,1) = correlation_matrix(1,1);
    a = linsolve(correlation_matrix,rho);    
end
