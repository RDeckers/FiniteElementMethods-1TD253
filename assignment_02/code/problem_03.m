function [ ] = problem_03( )
  [p,e,t] = create_mesh(0.1);
  [A,R, b, r] = assemble(p, e, t, 1e6, @f_sin, @g_const, @kappa);
  z = solve(A, R, b, r);
  for i = 1:size(p,2)
    x = p(1,i);
    y = p(2,i);
    z(i) = abs(z(i)-  exact(x,y));
  end
  figure(1)
  make_plot(p,e,t,z);
end

function z = exact(x, y)
  z = sin(pi*x)*sin(pi*y);
end