function [  ] = part_01(  )
  Hlog = [];
  Normlog = [];
  for i = 1:0.25:5
    h = 2^-i;
    [p,e,t] = create_mesh(h);
    [A,b] = assemble(p, e, t, @f);
    z = A\b;
    err = compute_delta(p, z);
    norm = sqrt(err.'*A*err);
    Hlog = [Hlog -i];
    Normlog = [Normlog log(norm)];
  end
  figure(1)
  plot(Hlog, Normlog)
  order = polyfit(Hlog, Normlog,1)
  hold on;
  fh = @(x) order(1)*x+order(2);
  ezplot(fh, [min(Hlog), max(Hlog)]);
  
  h = 1/8;
  [p,e,t] = create_mesh(h);
  [A,b] = assemble(p, e, t, @f);
  z = A\b;
  figure(2)
  make_plot(p,e,t,z)
  exact = compute_exact(p);
  figure(3)
  make_plot(p,e,t,exact);
  err = compute_delta(p, z);
  figure(4)
  make_plot(p,e,t,err);
  
  h = 1/16;
  [p,e,t] = create_mesh(h);
  [A,b] = assemble(p, e, t, @f);
  z = A\b;
  figure(5)
  make_plot(p,e,t,z)
  exact = compute_exact(p);
  figure(6)
  make_plot(p,e,t,exact);
  err = compute_delta(p, z);
  figure(7)
  make_plot(p,e,t,err);
end

function [ z ] = f( x, y )
  z=8*pi^2*sin(2*x*pi)*sin(2*y*pi);
end

function [ p,e,t ] = create_mesh( dx )
 R1 = [3,4, -1,1,1,-1, -1,-1,1,1]'; %unit square
 dl = decsg(R1); %make geometry
 [p, e, t] = initmesh(dl,'hmax', dx); %create a mesh with max segment size dx.
end

function [z] = exact(x, y)
  z = sin(2*pi*x)*sin(2*pi*y);
end

function e = compute_exact(points)
size(points)
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
