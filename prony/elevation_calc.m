function elevation = elevation_calc(signal,lambda,d)
    shape = size(signal);
    N = shape(1);
    phi = zeros(1,N-1);

    for i=2:N
        A1 = abs( signal(i-1) ) ;
        A2 = abs( signal(i)   );
        phi(i-1) = atan(  angle( A1/A2 * signal(i)/signal(i-1) )   );
    end

    elevation = zeros(1,N-1);
    for i=1:N-1
        elevation(i) = asin( lambda/d * phi(i)/(2*pi) );
    end
end

