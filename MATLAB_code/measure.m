function [y] = measure(x,pf)
% measure density
if size(x,3)>1 % if x is 3D matrix (vectorized with Np as the third dimension)
    y = x(:,pf.meas_pt,:);
else
    y = x(:,pf.meas_pt);
end
y(y<0) = 0;
end