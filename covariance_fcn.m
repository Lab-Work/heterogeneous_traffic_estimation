function R = covariance_fcn(tau,l)
% autocovariance function
% K(x,x') = exp(-|Tau|/l), where tau is the relative distance and l is the
% characteristic length scale
R = exp(-abs(tau)/l); % Ornstein?Uhlenbeck
end

