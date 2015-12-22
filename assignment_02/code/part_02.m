function [  ] = part_02(  )
  %try different h, 1/8, 1/16, 1/32? with different relative dt. plot norm.
  %problem_02_04()
%   h = 1/16;
%   dt = 0.5*h/(2*pi);
%   [p,e,t] = create_mesh(h);
%   size(p)
%   u = initialize(p);
%   [A, M, b] = assemble_2(p, e, t, @(x,y) 0, @beta_field, 0);
%   Fp = M+dt/2*A;
%   Fm = M-dt/2*A;
%   b = dt/2*b;
%   figure(1)
%   make_plot(p,e,t, u)
%   figure(2)
%   for i = dt:dt:1+dt/2
%     u = Fp\(Fm*u+b);
%     u(e(1,:)) = 0; %force the dirchlet boundary condition.
%     i
%   end
%   make_plot(p,e,t, u)
%   L2_norm(exact_error(p,u), M)
%   %F = pdeInterpolant(p,t,u);
  %interp = @(x,y) evaluate(F,[x;y]);
  %figure(3)
  %interp(2,3)
  %interp(1,-0)
  %ezcontourf(interp, [-1 1 -1 1])
  %axis equal
end

function [] = problem_02_04()
  H = [1/4 1/8 1/16 1/32];
  h_ref = 1/64;
  [p,e,t] = create_mesh(h_ref);
  u = initialize(p);
  [A, M, b] = assemble_2(p, e, t, @(x,y) 0, @beta_field, 0.1);
  dt = 0.5*h_ref/(2*pi);
  Fp = M+dt/2*A;
  Fm = M-dt/2*A;
  for i = dt:dt:1+dt/2
    u = Fp\(Fm*u+b);
    u(e(1,:)) = 0; %force the dirchlet boundary condition.
    i
  end
  Err = []
  F = pdeInterpolant(p,t,u);
  interp = @(x,y) evaluate(F,[x;y]);
  for h = H
    [p,e,t] = create_mesh(h_ref);
    u = initialize(p);
    [A, M, b] = assemble_2(p, e, t, @(x,y) 0, @beta_field, 0.1);
    dt = 0.5*h/(2*pi);
    Fp = M+dt/2*A;
    Fm = M-dt/2*A;
    
    for i = dt:dt:1+dt/2
      u = Fp\(Fm*u+b);
      u(e(1,:)) = 0; %force the dirchlet boundary condition.
      i
    end
    Err = [Err log(L2_norm(abs(approximate_error(p,u, interp)),A))]
  end
  figure(3)
  plot(log(H), Err)
  hold on;
  order = polyfit(log(H), Err,1)  % = 1.1894   -0.3067
  fh = @(x) x*order(1)+order(2);
  ezplot(fh, [min(log(H)), max(log(H))]);
end

function [ z ] = f( x, y )
  z=8*pi^2*sin(2*x*pi)*sin(2*y*pi);
end

function [ p,e,t ] = create_mesh( dx )
 g = @circleg;
 [p, e, t] = initmesh(g,'hmax', dx); %create a mesh with max segment size dx.
end

function L = L2_norm(e, A)
  L = sqrt(e'*A*e);
end
  
function [e] = approximate_error(points, solution, approximation)
  e = zeros(size(points,1),1);
  for i = 1:size(points,2)
      x = points(1,i);
      y = points(2,i);
      e(i) = approximation(x,y)-solution(i);
  end
end

function [e] = exact_error(points,solution)
  e = zeros(size(points,1),1);
  for i = 1:size(points,2)
      x = points(1,i);
      y = points(2,i);
      e(i) = u0_field(x,y)-solution(i);
  end
end

function u0 = initialize(points)
u0 = zeros(size(points,1),1);
  for i = 1:size(points,2)
      x = points(1,i);
      y = points(2,i);
      u0(i) = u0_field(x,y);
  end
end

function z = u0_field(x, y)
  x0 = 0.3;
  y0 = 0;
  r0 = 0.25;
  z = 1/2*(1-tanh(((x-x0)^2+(y-y0)^2)/r0^2-1));
end

function e = compute_exact(points)
  e = zeros(size(points,1),1);
  for i = 1:size(points,2)
      x = points(1,i);
      y = points(2,i);
      e(i) = exact(x,y);
  end
end

function [e] = compute_delta(points, solution)
  e = zeros(size(solution));
  for i = 1:size(points,2)
      x = points(1,i);
      y = points(2,i);
      e(i) = abs(exact(x,y)- solution(i));
  end
end

function [  ] = make_plot( p, e, t, z )
  pdeplot(p,e,t, 'xydata', z, 'mesh', 'on');
end

function [ x_out, y_out ] = beta_field( x, y )
 x_out = -2*pi*y;
 y_out = 2*pi*x;
end

