function my_first_fem_solver(N)
a = -1; % left end point of interval
b = 1; % right
h = (b-a)/N; % mesh size
x = a:h:b; % node coords
A=my_stiffness_matrix_assembler(x);
B=my_load_vector_assembler(x);
xi = A\B; % solve system of equations
%plot(x,xi) % plot solution
