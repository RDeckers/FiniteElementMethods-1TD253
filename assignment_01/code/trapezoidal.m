function [ rho ] = trapezoidal(x, delta)
    N = length(delta)+1;
    rho = zeros(N,1);
    rho(1) = (x(2)-x(1))^2*delta(1); %first interval, project backwards
    for i = 2:N-1
        h = x(i+1)-x(i);
        rho(i) = 0.5*h^2*(delta(i-1)+delta(i));
    end
    rho(N) = (x(N+1)-x(N))^2*delta(N-1); %last interval, project forwards
end

