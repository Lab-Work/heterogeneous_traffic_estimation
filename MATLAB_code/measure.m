function [y] = measure(x,pf)
% measure density
    y = x(:,pf.meas_pt);
    y(y<0) = 0;
end