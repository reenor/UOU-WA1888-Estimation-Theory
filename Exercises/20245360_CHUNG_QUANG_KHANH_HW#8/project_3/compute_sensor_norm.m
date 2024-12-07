function [ ynorm ] = compute_sensor_norm( y , mode )
% function [ ynorm ] = compute_sensor_norm( y , mode )
% mode == 1 : norm
% mode == 2 : difference norm

if (nargin < 2 )
  mode = 1;
end


N = size(y,2);
ynorm = zeros(1,N);


if ( mode == 1 )
  for i = 1:N
    ynorm(i) = norm( y(:,i) );
  end
elseif ( mode == 2 )
  for i = 2:N
    ynorm(i) = norm( y(:,i) - y(:,i-1));
  end
else
  for i = 2:N
     ynorm(i) = norm(y(:,i)) - norm(y(:,i-1));
  end
end