function [y] = measure_true(x,pf)
% measure density
    y = x(:,pf.meas_pt);
    y = y+normrnd(0,0.03*ones(2,length(pf.meas_pt)));
    y(y<0) = 0; % non-negativity
end