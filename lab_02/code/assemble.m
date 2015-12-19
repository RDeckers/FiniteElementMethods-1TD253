function [A,R,b,r] = assemble( p,e,t, gamma, f, g, kappa)

N = size(p,2);
A = sparse(N,N);
R = sparse(N,N);
b = zeros(N,1);
r = zeros(N,1);

for K = 1:size(t,2); % loop over the triangles
  nodes = t(1:3,K); % find triangle K's nodes
  
  x = p(1,nodes);
  y = p(2,nodes);
  
  [AK, bK] = create_AK_bK(x,y, f, kappa);
  
  A(nodes,nodes) = A(nodes,nodes)+AK; % add AK(i,j), i,j=1,2,3, to A(nodes(i),nodes(j))
  b(nodes) = b(nodes) + bK;
end
for E = 1:size(e,2)
  nodes = e(1:2,E);
  x = p(1,nodes);
  y = p(2,nodes);
  length_E = sqrt((x(1)-x(2))^2+(y(1)-y(2))^2);
  R(nodes,nodes) = R(nodes,nodes) + gamma*length_E*[2 1;1 2]/6;
  r(nodes) = r(nodes) + gamma*length_E*g(x,y)'/2;
end

end

