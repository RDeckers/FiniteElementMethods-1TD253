function [  ] = make_plot( p, e, t, z )
  pdeplot(p,e,t, 'xydata', z, 'mesh', 'on');
end

