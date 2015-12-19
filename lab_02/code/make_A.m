function [ A ] = make_A( dx )
%MAKE_A Summary of this function goes here
%   Detailed explanation goes here
R1 = [3,4, 0,1,1,0, 0,0,1,1]'; %unit square
[dl, bt] = decsg(R1); %make geometry
[p, e, t] = initmesh(dl,'Hmax', dx) %create a mesh with max segment size dx.
for K = 1:size(t,2); % loop over the triangles
  nodes = t(1:3,K); % find triangle K's nodes
  
  x = p(1,nodes);
  y = p(2,nodes);
  area_K = polyarea(x,y);
  
  AK =  [ones(3); x; y]% compute the (3 x 3) stiffness matrix AK
  A(nodes,nodes) = A(nodes,nodes)+AK; % add AK(i,j), i,j=1,2,3,
% to A(nodes(i),nodes(j))
end
end

