function my_first_fem_solver(N)
a = 0; % left end point of interval
b = 1; % right
h = (b-a)/N; % mesh size
x = a:h:b; % node coords
low =  0;
high = 0.;
B=my_load_vector_assembler(x, low, high, @f);
M = mass_matrix(x); %last points kept fixed
%M(1,:) = [];
%M(end,:) = [];
A_fixed =stiffness_matrix_fixed(x,low, high);
xi_fixed = A_fixed\B; % solve system of equations
reversed = (-A_fixed*xi_fixed);
lap_fixed = M \(reversed(2:end-1));
figure(1)
plot(x,xi_fixed) % plot solution
figure(2)
hold on;
F = arrayfun(@f, x(2:end-1)).';
delta = abs(F+lap_fixed);
plot(x(2:end-1), delta)
x
rho = trapezoidal(x, abs(F+lap_fixed));
figure(3)
hold on;
plot(rho)
%figure(3)
%plot(x, lap_fixed-lap_huh)

function y=f(x)
  %y=2;
  %y=exp(-1000*(x-0.5)^2);
  %y = pi^2*49*sin(x*pi*7);
  y = x*(x-1)