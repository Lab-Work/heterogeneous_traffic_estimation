function K = covariance(pf,den)
% covariance matrix K

% tau: distance between two measurements
% tau = (pf.meas_pt(end)-pf.meas_pt(1))/(length(pf.meas_pt)-1);
% n = 2*length(pf.meas_pt); % number of measurements
% inverse_eye = (ones(n)-eye(n));
% inverse_eye(inverse_eye==1)=Inf;
% tau_m = tau*inverse_eye;
% 
% K = pf.meas_stdev^2 .* covariance_fcn(tau_m,den);

% autocovariance matrix R
tau = (pf.meas_pt(end)-pf.meas_pt(1))/(length(pf.meas_pt)-1);
tau_m = tau*[0 inf inf inf inf inf
    inf 0 inf inf inf inf
    inf inf 0 inf inf inf
    inf inf inf 0 inf inf
    inf inf inf inf 0 inf
    inf inf inf inf inf 0];

K = pf.meas_stdev^2 .* covariance_fcn(tau_m,den);

end

