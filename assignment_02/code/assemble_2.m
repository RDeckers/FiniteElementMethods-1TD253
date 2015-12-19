function [Fp, Fm,b] = assemble_2( p,e,t, f, beta, epsilon, delta_t)

I = eye(length(p));
N = size(p,2);
Fp = sparse(N,N);
Fm = sparse(N,N);
b = zeros(N,1);

for K = 1:size(t,2); % loop over the triangles
  nodes = t(1:3,K); % find triangle K's nodes
  
  x = p(1,nodes);
  y = p(2,nodes);
  
  [AK, bK, BK] = create_AK_bK2(x,y, epsilon, f, beta);
  Fp(nodes,nodes) = Fp(nodes,nodes) + BK+delta_t/2*AK;
  Fm(nodes,nodes) = Fm(nodes,nodes) + BK-delta_t/2*AK;
  % add AK(i,j), i,j=1,2,3, to A(nodes(i),nodes(j))
  b(nodes) = b(nodes) + delta_t/2*bK;
  
end
Fp(e(1,:),:) = I(e(1,:),:); % replace the rows corresponding
Fm(e(1,:),:) = I(e(1,:),:); % replace the rows corresponding
% to the boundary nodes by corresponding
% rows of I
b(e(1,:)) = 0; % put the boundary value into the RHS
end




