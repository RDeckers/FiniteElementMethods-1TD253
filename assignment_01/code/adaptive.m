function adaptive(N, N_max, lambda)
a = -1; % left end point of interval
b = 1; % right
h = (b-a)/N; % mesh size
x = a:h:b; % node coords
low =  0;
high = 0.;
figure(4)
hold on;
figure(3)
hold on;
figure(2)
hold on;
figure(1)
hold on;
N_array = [];
total_residual_array =  [];

while N < N_max
    N_array = [N_array N];
    B=my_load_vector_assembler(x, low, high, @f);
    M = mass_matrix(x); %last points kept fixed
    A_fixed =stiffness_matrix_fixed(x,low, high);
    xi_fixed = A_fixed\B; % solve system of equations
    reversed = (-A_fixed*xi_fixed);
    lap_fixed = M \(reversed(2:end-1));
    F = arrayfun(@f, x(2:end-1)).';
    delta = abs(F+lap_fixed);
    rho = trapezoidal(x, delta);
    total_residual_array = [total_residual_array sqrt(sum(rho.*rho))];
    threshold = max(rho)*lambda;
    for i = 1:length(rho)
        if rho(i) > threshold
            x = [x (x(i+1)+x(i))/2];
        end
    end
    x = sort(x);
    N = size(x,2)
end
%total_residual_array(end)
%rho
N_array = [N_array N]
B=my_load_vector_assembler(x, low, high, @f);
M = mass_matrix(x); %last points kept fixed
A_fixed =stiffness_matrix_fixed(x,low, high);
xi_fixed = A_fixed\B; % solve system of equations
reversed = (-A_fixed*xi_fixed);
lap_fixed = M \(reversed(2:end-1));
F = arrayfun(@f, x(2:end-1)).';
delta = abs(F+lap_fixed);
rho = trapezoidal(x, delta);
total_residual_array = [total_residual_array sqrt(sum(rho.*rho))];
figure(1)
plot(x, xi_fixed)
figure(2)
semilogy(x(2:end), [1./diff(x)])
figure(3)
loglog(N_array, total_residual_array)
figure(4)
plot(x(1:end-1), rho)
%x
function y=f(x)
%y=2;
y=exp(-1000*x^2)+10^-3;
%y = pi^2*49*sin(x*pi*7);
%y = x*(x-1);