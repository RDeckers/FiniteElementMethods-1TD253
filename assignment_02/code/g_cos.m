function [ z ] = g_cos( x, y )
  if(x == 0)
    z = cos(2*pi*y);
  else
    z = 0;
  end
end

