function [ z ] = solve( A, R, b, r )
%SOLVE Summary of this function goes here
%   Detailed explanation goes here
  z = (A+R)\(b+r);
  %pdeplot(p,e,t,'xydata',z, 'mesh', 'on');
end

