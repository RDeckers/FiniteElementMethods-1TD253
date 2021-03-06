function A=stiffness_matrix_fixed(x, low, high)
%
% Returns the assembled stiffness matrix A.
% Input is a vector x of node coords.
%
N = length(x) - 1;
% number of elements
A = zeros(N+1, N+1); % initialize stiffnes matrix to zero
for i = 1:N
% loop over elements
h = x(i+1) - x(i); % element length
n = [i i+1];
% nodes
A(n,n) = A(n,n) + [1 -1; -1 1]/h; % assemble element stiffness
end
A(1,:) = zeros(N+1,1);
A(1,1) = 1; % adjust for BC
A(N+1, :) = zeros(N+1,1);
A(N+1,N+1) = 1; % adjust for BC
%A=sparse(A); %Vrrrrrrooooom!