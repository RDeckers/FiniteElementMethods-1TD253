function [A, M, b] = assemble_2( p,e,t, f, beta, epsilon)

I = eye(length(p));
N = size(p,2);
A = sparse(N,N);
M = sparse(N,N);
b = zeros(N,1);

for K = 1:size(t,2); % loop over the triangles
  nodes = t(1:3,K); % find triangle K's nodes
  
  x = p(1,nodes);
  y = p(2,nodes);
  
  [AK, bK, MK] = create_AK_bK2(x,y, epsilon, f, beta);
  A(nodes,nodes) = A(nodes,nodes) + AK;
  M(nodes,nodes) = M(nodes,nodes) + MK;
  % add AK(i,j), i,j=1,2,3, to A(nodes(i),nodes(j))
  b(nodes) = b(nodes)+ bK;
end
%A(e(1,:),:) = I(e(1,:),:); % replace the rows corresponding
%M(e(1,:),:) = I(e(1,:),:); % replace the rows corresponding
% to the boundary nodes by corresponding
% rows of I
b(e(1,:)) = 0; % put the boundary value into the RHS
end




