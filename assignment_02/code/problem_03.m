function [order] = problem_03()
  X =[];
  Norm = [];
  %loop over different grid sizes
  %in logspace
  for j = 0.5:0.1:6.5 
    dx = 2^-j;
    [p,e,t] = create_mesh(dx);
    [A,R, b, r] = assemble(p, e, t, 1e6, @f_sin, @g_const, @kappa);
    z = solve(A, R, b, r);
    for i = 1:size(p,2)
      x = p(1,i);
      y = p(2,i);
      z(i) = abs(z(i)-  exact(x,y));
    end
    norm = log(z.'*A*z);
    Norm = [Norm norm];
    X = [X log(dx)];
  end
  figure(1)
  plot(X, Norm)
  order = polyfit(X, Norm,1);
  hold on;
  fh = @(x) order(1)*x+order(2);
  ezplot(fh, [-7, 0]);
  
  dx = 0.05;
  [p,e,t] = create_mesh(dx);
  [A,R, b, r] = assemble(p, e, t, 1e6, @f_sin, @g_const, @kappa);
  z = solve(A, R, b, r);
  figure(2)
  make_plot(p,e,t, z);
  
end


function z = exact(x, y)
  z = sin(pi*x)*sin(pi*y);
end