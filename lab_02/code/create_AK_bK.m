function [ AK, bK ] = create_AK_bK( x, y, f, kappa)
  %x, y are triplets of vertices, f and kappa are funciton handles.
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
  %and compute kappa and f at the centroid
  AK = (b.'*b+c.'*c)*kappa(x_c, y_c)*area_K;
  bK = f(x_c, y_c)*area_K/3;
end
