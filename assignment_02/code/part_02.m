function [  ] = part_02(  )
  h = 1/16;
  dt = 0.5*h/(2*pi);
  [p,e,t] = create_mesh(h);
  size(p)
  u = initialize(p);
  [Fp, Fm, b] = assemble_2(p, e, t, @(x,y) 0, @beta_field, 0, dt);

  figure(1)
  make_plot(p,e,t, u)
  figure(2)
  for i = dt:dt:1+dt/2
    u = Fp\(Fm*u+b);
    i
  end
  make_plot(p,e,t, u)
end

function [ z ] = f( x, y )
  z=8*pi^2*sin(2*x*pi)*sin(2*y*pi);
end

function [ p,e,t ] = create_mesh( dx )
 g = @circleg;
 [p, e, t] = initmesh(g,'hmax', dx); %create a mesh with max segment size dx.
end

function [z] = exact(x, y)
  z = sin(2*pi*x)*sin(2*pi*y);
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
  pdeplot(p,e,t, 'xydata', z, 'zdata', z, 'mesh', 'on');
end

function [ x_out, y_out ] = beta_field( x, y )
 x_out = -2*pi*y;
 y_out = 2*pi*x;
end

