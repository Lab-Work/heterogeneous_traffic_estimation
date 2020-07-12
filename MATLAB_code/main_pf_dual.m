close all
clear all
% October 2019
% Yanbing Wang
% Vanderbilt University

% parameter-adaptive particle filter for heterogeneous traffic state
% estimation (vectorized implementation to speed up runtime)


%% =========== test setup ==============
% model parameters
model_true = model_param; model_true.model_name = 'True model';
model_est = model_param; model_est.model_name = 'Creeping model';

% test parameters
perturbation = [1]; % for sensitivity analysis
test = 1; % 1: overtaking, 2: queue clearance, 3: congested flow, 4:creeping
% Npc = [500 800 1000 1500 2000]; % number of particles
Npc = 500;
Nr = 1; % number of simulation runs
spatial_correlation = true; % toggle to turn on/off spatial correlation
show_sim = false; % toggle to turn on/off simulation display
show_est = false; % toggle to turn on/off estimation display
len = 60; % characteristic length scale for spatial correlation
T_pi = []; % particle influence summary table

% save data
directory = pwd; % create a new folder to save test run
foldername = sprintf('test%d_fil_papf_corr',test);
mkdir(foldername);
directory = fullfile(directory,foldername);

% PF parameters
pf = params_PF;

%% initial conditions and boundary conditions
for pc = 1:length(Npc)
    % for pc = 1:length(perturbation)
    pf.Np = Npc(pc);
    fprintf('#particles %d\n',pf.Np);
    for r = 1:Nr
        
        fprintf('run %d\n',r);
        
        % initialize states based on test setup
        U0_true = initialize(model_true, test);
        U0_est = initialize(model_est, test);
        U_true_c1 = zeros([model_true.M model_true.N]); U_true_c2 = U_true_c1;
        U_est_c1 = U_true_c1; U_est_c2 = U_true_c1; U_sim_c1 = U_true_c1; U_sim_c2 = U_true_c1;
        U_true_c1(1,:) = U0_true(1,:); U_true_c2(1,:) = U0_true(2,:);
        U_est_c1(1,:) = U0_est(1,:); U_est_c2(1,:) = U0_est(2,:);
        
        % simulation of the true state and the estimator
        pf.meas_pt = [3 model_est.N/2 model_est.N-3];
        pf.init_stdev = 0.05; % 0.1
        pf.model_stdev = [0.02,0.03];% test1 0.02,0.03 test2: 0.05,0.05
        pf.meas_stdev = 0.06; % 0.06
        
        U_true = cell([1 model_true.M]); U_est = cell([1 model_true.M]);
        U_meas_true = cell([1 model_true.M]); U_res = cell([1 model_true.M]);
        U_true{1} = U0_true; U_est{1} = U0_est;
        U_meas_true{1} = measure_true(U_true{1},pf);
        
        for n = 1:model_true.M-1
            theta_true = [model_true.vm1, model_true.rm1, model_true.rm2];
            U_true{n+1} = solver_flow_dual(n,U_true{n},model_true,theta_true);
            U_true_c1(n+1,:) = U_true{n+1}(1,:);
            U_true_c2(n+1,:) = U_true{n+1}(2,:);
            
            U_est{n+1} = solver_flow(n,U_est{n},model_est);
            U_sim_c1(n+1,:) = U_est{n+1}(1,:);
            U_sim_c2(n+1,:) = U_est{n+1}(2,:);
            U_meas_true{n+1} = measure_true(U_true{n+1},pf);
            
            % **************** plot forward sim *****************
            if show_sim
                if mod(n,5)==0
                    fig = plot_compare(n,U_true,U_est,model_est);
                    %             filename = sprintf('simulation_%d_%03d',test,n);
                    %             path = fullfile(directory,filename);
                    %             saveas(gca,path,'png')
                end
                drawnow
            end
        end
        
        %% initialize particles and weights
        R = covariance(pf,100); % measurement covariance matrix
        x = zeros([size(U0_est),pf.Np]); % dimension: classes x cells x particles
        wt = ones(pf.Np, 1)/pf.Np; % initialize particles to have uniform weights
        y_miu = zeros(size(measure(U_est{1},pf))); % measurement error mean
        m = numel(y_miu); % measurement vector size
        magic1 = (2 * pi)^(-m/2) * (sqrt(sum(sum(abs(R).^2))))^(-1/2); % for computing particle weights
        magic2 = inv(R)*(-0.5); % for computing particle weights
        sum_init_1 = 0; sum_init_2 = 0;
        
        % ******************* spatial correlation **********************
        if spatial_correlation % correlation parameters
            tau = 0:1:model_est.N-1;
            tau_m = triu(toeplitz(tau));
            tau_m = tau_m-tau_m'; % distance in # cells
            R_m = covariance_fcn(tau_m,len);
            [V,D] = eig(R_m);
            x_miu = zeros(size(tau));
            init_miu = x_miu;
            
            %^^^^^^^^^^^^^^ initial noise vectorized ^^^^^^^^^^^^^^^^^^^^^^^
            for i = 1:size(V,1)
                sum_init_1 = sum_init_1 + sqrt(D(i,i)) * randn(pf.Np,1) * pf.init_stdev * V(:,i)';
                sum_init_2 = sum_init_2 + sqrt(D(i,i)) * randn(pf.Np,1) * pf.init_stdev * V(:,i)';
            end
            x = cat(3, (U0_est(1,:) + sum_init_1)', (U0_est(2,:) + sum_init_2)');
            x = permute(x, [3 1 2]);
            x(x<0)=0;
            
            % ******************* Uncorrelated vectorized **********************
        else
            x = U0_est + normrnd(0,pf.init_stdev,size(x)); % state-independent noises
        end
        
        %% Perform time propagation to obtain priori particles
        tic
        x_next = x;
        flow_next_param = x;
        y_next = zeros(size(measure(U_est{1},pf)));
        N_eff = zeros(1,model_true.M);
        U_res{1} = U0_est;
        theta = [model_est.vm1, model_est.rm1, model_est.rm2]'; % theta^{0|0}
        Nm = pf.Np;
        theta_matrix = zeros(3,Nm);
        
        %   ****************************** PARTICLE FILTER STARTS ***********************************
        for n = 2:model_true.M
            %     fprintf('running time = %d\n',n);
            U_est_temp = mean(x,3);
            x(x<0)=0;
            y_true = [U_meas_true{n}(1,:)';U_meas_true{n}(2,:)']; % true measurement
            
            % boundaries
            d1l_temp = model_est.d1l(n);
            d2l_temp = model_est.d2l(n);
            d1r_temp = model_est.d1r(n);
            d2r_temp = model_est.d2r(n);
            
            % =============== dual parameter estimation ========================
            
            switch test
                case 1
                    theta_stdev = [0.005, 0.005, 0.005];
                case 2
                    theta_stdev = [0.01, 0.01, 0.005];
                case 3
                    theta_stdev = [0.005, 0.005, 0.005];
                case 4
                    theta_stdev = [0.005, 0.002, 0.002];
            end
            
            
            % ^^^^^^^^^^^^^^ process noise vectorized ^^^^^^^^^^^^^^^^^^^^^^^^^^^
            if spatial_correlation
                noise_x_1 = 0; noise_x_2 = 0; sum_y = 0;
                for i = 1:size(V,1)
                    noise_x_1 = noise_x_1 + sqrt(D(i,i)) * randn(pf.Np,1) * pf.model_stdev(1) * V(:,i)';
                    noise_x_2 = noise_x_2 + sqrt(D(i,i)) * randn(pf.Np,1) * pf.model_stdev(2) * V(:,i)';
                end
                noise_x = cat(3, noise_x_1, noise_x_2);
                noise_x = permute(noise_x, [3 2 1]);
            else
                noise_x = normrnd(0,pf.model_stdev(1),size(x_next));
            end
            
            %------ predict vectorized ---------
            theta_matrix = repmat(theta(:,n-1),1,Nm) + normrnd(zeros([size(theta,1),Nm]),repmat(theta_stdev',1,Nm)); % theta^{k|k-1}
            
            for p = 1:Nm
                %------ predict ---------
                flow_next_param(:,:,p) = solver_flow_dual(n-1,U_est_temp,model_est,theta_matrix(:,p)); % x^{k|k-1}
                flow_next_param(flow_next_param<0) = 0;
                
            end
            
            %------ predict noise vectorized ---------
            flow_next_param = flow_next_param + noise_x;
            
            %------- parameter update vectorized -----
            noise_y = y_miu + normrnd(0,pf.meas_stdev,[size(y_miu), Nm]);
            y_next = measure(flow_next_param,pf) + noise_y;
            y_next(y_next<0) = 0;
            y_n = reshape(permute(y_next,[2 1 3]), m, Nm);
            
            wtp = magic1 * ...
                exp(sum((y_true - y_n).* (magic2 * (y_true - y_n))));
            
            %------ normalize weight vectorized ----
            wtp = wtp./sum(wtp);
            N_eff_param(n) = (sum(wtp.^2))^(-1);
            
            
            % ------ resampling -----
            for p = 1:Nm
                theta_matrix(:,p) = theta_matrix(:,(find(rand <= cumsum(wtp),1)));
            end
            
            theta(:,n) = mean(theta_matrix,2); % theta^{k|k}
            
            
            
            % =============== dual state estimation ========================
            
            
            %^^^^^^^^^^^^^^ process noise vectorized ^^^^^^^^^^^^^^^^^^^^^^^^^^^
            if spatial_correlation
                noise_x_1 = 0; noise_x_2 = 0; sum_y = 0;
                for i = 1:size(V,1)
                    noise_x_1 = noise_x_1 + sqrt(D(i,i)) * randn(pf.Np,1) * pf.model_stdev(1) * V(:,i)';
                    noise_x_2 = noise_x_2 + sqrt(D(i,i)) * randn(pf.Np,1) * pf.model_stdev(2) * V(:,i)';
                end
                noise_x = cat(3, noise_x_1, noise_x_2);
                noise_x = permute(noise_x, [3 2 1]);
            else
                noise_x = normrnd(0,pf.model_stdev(1),size(x_next));
            end
            
            %^^^^^^^^^^^^^^ measurement noise vectorized ^^^^^^^^^^^^^^^^^^^^^^^
            noise_y = y_miu + normrnd(0,pf.meas_stdev,[size(y_miu), pf.Np]);
            
            for p = 1:pf.Np
                
                %^^^^^^^^^^^^^^ boundary noise ^^^^^^^^^^^^^^^^^^^^^^^^^^
                if test == 2
                    model_est.d1l(n) = d1l_temp + wblrnd(0.5,6);
                    model_est.d2l(n) = d2l_temp + wblrnd(0.3,4);
                    model_est.d1r(n) = d1r_temp + wblrnd(0.06,4);
                    model_est.d2r(n) = d2r_temp + wblrnd(0.06,4);
                else
                    model_est.d1l(n) = d1l_temp + wblrnd(0.06,4);
                    model_est.d2l(n) = d2l_temp + wblrnd(0.06,4);
                    model_est.d1r(n) = d1r_temp + wblrnd(0.06,4);
                    model_est.d2r(n) = d2r_temp + wblrnd(0.06,4);
                end
                
                %------ model predict --------------
                flow_next(:,:,p) = solver_flow_dual(n-1,x(:,:,p),model_est,theta(:,n));
                
            end
            %------ predict noise vectorized ---------
            x_next = flow_next + noise_x;
            x_next(x_next<0)=0;
            
            %------- parameter update vectorized -----
            %             noise_y = y_miu + normrnd(0,pf.meas_stdev,[size(y_miu), Nm]);
            y_next = measure(flow_next,pf) + noise_y;
            y_next(y_next<0) = 0;
            y_n = reshape(permute(y_next,[2 1 3]), m, pf.Np);
            
            wt = magic1 * ...
                exp(sum((y_true - y_n).* (magic2 * (y_true - y_n))));
            
            %------ normalize weight vectorized ----
            wt = wt./sum(wt);
            N_eff(n) = (sum(wt.^2))^(-1);
            
            
            % ------ resampling -----
            for p = 1:pf.Np
                x_next(:,:,p) = x_next(:,:,(find(rand <= cumsum(wt),1)));
            end
            
            
            U_est_temp = mean(x_next,3); % take mean
            %     U_est_temp = x_next(:,:,ind); % take MAP
            U_est_temp(U_est_temp<0) = 0;
            U_res{n} = U_est_temp;
            U_est_c1(n,:) = U_est_temp(1,:);
            U_est_c2(n,:) = U_est_temp(2,:);
            
            
            %     --------------------------- plot and save ---------------------------
            if show_est
                if mod(n,5)==0
                    fig = plot_est(n,U_true,U_res,model_est,pf,U_meas_true,x_next);
                    drawnow
                    %                     filename = sprintf('pf_test%d_%03d',test,n);
                    %                     path = fullfile(directory,filename);
                    %                     saveas(gca,path,'png')
                end
            end
            
            x = x_next;
        end %n
        runtime(r) = toc;
        
        %% write video
        % write_video(directory,foldername);
        
        %% plot estimated densities
        % plot_six(U_true_c1,U_sim_c1,U_est_c1,U_true_c2,U_sim_c2,U_est_c2,model_true);
        % plot_six_gray(U_true_c1,U_sim_c1,U_est_c1,U_true_c2,U_sim_c2,U_est_c2,model_true);
        
        %% plot effective particle size
        % plot_eps(N_eff)
        N_eff_mean(r) = mean(N_eff);
        
        %% compute MAE
        % between simulated and true
        MAE_sim_c1 = sum(sum(abs(U_true_c1 - U_sim_c1)))/(model_true.N*model_true.M);
        MAE_sim_c2 = sum(sum(abs(U_true_c2 - U_sim_c2)))/(model_true.N*model_true.M);
        % between estimated and true
        MAE_est_c1 = sum(sum(abs(U_true_c1 - U_est_c1)))/(model_true.N*model_true.M);
        MAE_est_c2 = sum(sum(abs(U_true_c2 - U_est_c2)))/(model_true.N*model_true.M);
        MAE(:,r) = [MAE_sim_c1;MAE_sim_c2;MAE_est_c1;MAE_est_c2;(MAE_sim_c1-MAE_est_c1)/MAE_sim_c1*100;(MAE_sim_c2-MAE_est_c2)/MAE_sim_c2*100];
        if r == Nr
            disp(MAE);
        end
        
        %% compute BEEQ Bayesian estimation error quotient (BEEQ).
        beta_c1(r) = norm(U_true_c1-U_est_c1)/norm(U_true_c1-U_sim_c1);
        beta_c2(r) = norm(U_true_c2-U_est_c2)/norm(U_true_c2-U_sim_c2);
        
        % save(sprintf('vm%.2f_rm1%.2f_rm2%.2f',model_est.vm1,model_est.rm1,model_est.rm2));
        % save(sprintf('test%d_run%d',test,r));
        
    end % end of all runs
    log_beta_c1 = sum(log10(beta_c1))/Nr;
    log_beta_c2 = sum(log10(beta_c2))/Nr;
    true_beta_c1 = 10^log_beta_c1;
    true_beta_c2 = 10^log_beta_c2;
    
    %% save results to a table - influence of particles
    T_pi = [T_pi; Npc(pc) mean(runtime) true_beta_c1 true_beta_c2 mean(N_eff_mean)...
        mean(MAE(5,:)) mean(MAE(6,:))];
    % clear MAE runtime N_eff_mean
    
end % end of all particle choices
colNames = {'particles','runtime','beta_1','beta_2','mean_N_eff','improvement_1','improvement_2'};
T_pi = array2table(T_pi,'VariableNames',colNames);
save(sprintf('test%d_corr_adaptive_mac.mat',test),'T_pi');
T_pi