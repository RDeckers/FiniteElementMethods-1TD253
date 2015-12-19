function [ AK, bK, BK ] = create_AK_bK2( x, y, epsilon, f, beta)
  %x, y are triplets of vertices, f is a function handle.
  area_K = polyarea(x,y);
  %abc matrix
  Z = [ones(1,3); x; y].'; 
  %solve for the three abc vectors
  abc = [Z\[1;0;0] Z\[0;1;0] Z\[0;0;1]]; 
  b = abc(2,:); 
  c = abc(3,:); 
  %take the centroid coordinates
  x_c = mean(x);
  y_c = mean(y);
  %evaluate the given expression for AK and bK
  %compute f at the centroid
  BK = 1/12*area_K*[2 1 1; 1 2 1; 1 1 2];
  bK = f(x_c, y_c)*area_K/3;
  [beta_1, beta_2] = beta(x_c, y_c);
  beta_int = (beta_1*b+beta_2*c)*area_K/3; %same for all i
  AK = epsilon*(b.'*b+c.'*c)*area_K;
  AK = AK + [beta_int; beta_int; beta_int];
end
