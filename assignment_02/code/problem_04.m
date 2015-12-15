function [  ] = problem_04(  )
  dx = 0.05;
  [p,e,t] = create_mesh(dx);
  [A,R, b, r] = assemble(p, e, t, 1e6, @f_const, @g_cos, @kappa);
  z = solve(A, R, b, r);
  figure(3)
  make_plot(p,e,t, z);
end

