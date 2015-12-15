function [ AK, bK ] = create_AK_bK( x, y, f, kappa)
  area_K = polyarea(x,y);
  Z = [ones(1,3); x; y].'; %abc matrix
  abc = [Z\[1;0;0] Z\[0;1;0] Z\[0;0;1]];
  b = abc(2,:);
  c = abc(3,:);
  x_c = mean(x);
  y_c = mean(y);
  AK = (b.'*b+c.'*c)*kappa(x_c, y_c)*area_K;
  bK = f(x_c, y_c)*area_K/3;
end
