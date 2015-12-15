function [  ] = make_plot( p, e, t, z )
  pdeplot(p,e,t,'zdata',z, 'mesh', 'on');
end

