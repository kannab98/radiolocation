function R_xx = corr_mat(signal,p)
    R_xx = zeros(p,p);
    for i=1:p
        for j=1:p
            R_xx(i,j)= corr(signal,i-j);
        end
    end
end

