function M = mass_matrix( x )
    N = length(x) - 2;
    diag = zeros(N,1);
    h = zeros(N,1);
    upper = zeros(N-1,1);
    for i = 1:N+1
        h(i) = x(i+1)-x(i);
    end
    
    %diag(1) = (1/3*(x(2)^3-x(1)^3)+x(2)^2*h(1)-x(2)*(x(2)^2-x(1)^2))/h(1)^2;
    for i = 2:N+1
        diag(i-1) =           (1/3*(x(i+1)^3-x(i)^3)+x(i+1)^2*h(i)-x(i+1)*(x(i+1)^2-x(i)^2))/h(i)^2;
        diag(i-1) = diag(i-1) + (1/3*(x(i)^3-x(i-1)^3)+x(i-1)^2*h(i-1)-x(i-1)*(x(i)^2-x(i-1)^2))/h(i-1)^2;
    end
    %diag(N+1) = (1/3*(x(N+1)^3-x(N)^3)+x(N)^2*h(N)-x(N)*(x(N+1)^2-x(N)^2))/h(N)^2;
    
    for i = 1:N-1
        upper(i) = (-1/3*(x(i+1)^3-x(i)^3)+1/2*(x(i)+x(i+1))*(x(i+1)^2-x(i)^2)-x(i)*x(i+1)*h(i))/h(i)^2;
        %lower(i) = (-1/3*(x(i+1)^3-x(i)^3)+1/2*(x(i)+x(i+1))*(x(i+1)^2-x(i)^2)-x(i)*x(i+1)*h(i))/h(i)^2;
    end
    M = gallery('tridiag', upper, diag, upper);
end

