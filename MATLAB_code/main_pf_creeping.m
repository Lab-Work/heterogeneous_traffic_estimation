close all
% clc
% clear all
% October 2019
% Yanbing Wang
% Vanderbilt University

% =========== test setup ==============
test = 1; % 1: overtaking, 2: congested traffic, 3: queue clearance, 4:creeping traffic
Npc = 1500; % number of particles
Nr = 1; % number of simulation runs
spatial_correlation = true; % toggle SCNM
show_sim = false; % toggle simulation plots
show_est = false; % toggle estimation plots
len = 60; % 60 characteristic length for spatial correlation
T_pi = []; % particle influence summary table
directory = pwd;
foldername = sprintf('test%d_fil',test);
mkdir(foldername);
directory = fullfile(directory,foldername);

% =========== model parameters ==============
model_true = model_param; model_true.model_name = 'True model';
model_est = model_param; model_est.model_name = 'Creeping model';

% =========== PF parameters ==============
pf = params_PF;

for pc = 1:length(Npc)
    pf.Np = Npc(pc);
    fprintf('#particles %d\n',pf.Np);
    for r = 1:Nr
        fprintf('run %d\n',r);
        
% ==== initial conditions and boundary conditions ===
        U0_true = initialize(model_true, test);
        U0_est = initialize(model_est, test);  
        pf.meas_pt = [3 model_est.N/2 model_est.N-3]; % sensor location
        pf.init_stdev = 0.05; % initial noise stdev
        pf.model_stdev = [0.02,0.03]; % model noise stdev for class 1 and 2
        pf.meas_stdev = 0.06; % 0.06measurement noise stdev

        U_true_c1 = zeros([model_true.M model_true.N]); U_true_c2 = U_true_c1;
        U_est_c1 = U_true_c1; U_est_c2 = U_true_c1; U_sim_c1 = U_true_c1; U_sim_c2 = U_true_c1;
        U_true_c1(1,:) = U0_true(1,:); U_true_c2(1,:) = U0_true(2,:);
        U_est_c1(1,:) = U0_est(1,:); U_est_c2(1,:) = U0_est(2,:);
        
% ==== simulate true state and the predicted state ====
        U_true = cell([1 model_true.M]); U_est = cell([1 model_true.M]);
        U_meas_true = cell([1 model_true.M]); U_res = cell([1 model_true.M]);
        U_true{1} = U0_true; U_est{1} = U0_est;
        U_meas_true{1} = measure_true(U_true{1},pf);
        
        for n = 1:model_true.M-1
            U_true{n+1} = solver_flow(n,U_true{n},model_true);
            U_true_c1(n+1,:) = U_true{n+1}(1,:);
            U_true_c2(n+1,:) = U_true{n+1}(2,:);
            
            U_est{n+1} = solver_flow(n,U_est{n},model_est);
            U_sim_c1(n+1,:) = U_est{n+1}(1,:);
            U_sim_c2(n+1,:) = U_est{n+1}(2,:);
            U_meas_true{n+1} = measure_true(U_true{n+1},pf);
            
            % **************** plot forward sim *****************
            if show_sim
%                 if mod(n,5)==0
if n==5
                    fig = plot_compare(n,U_true,U_est,model_est);
                    filename = sprintf('simulation_%d_%03d',test,n);
                    path = fullfile(directory,filename);
                    saveas(gca,path,'png')
                end
                drawnow
            end
        end
% ================ end of simulation ================    
       
% ================ initialize particles ================
        R = covariance(pf,100); % measurement covariance matrix
        x = zeros([size(U0_est),pf.Np]); % class x cell x particles
        wt = ones(pf.Np, 1)/pf.Np;
        y_miu = zeros(size(measure(U_est{1},pf)));
        
        % ******************* spatial correlation **********************
        if spatial_correlation
            tau = 0:1:model_est.N-1;
            tau_m = triu(toeplitz(tau));
            tau_m = tau_m-tau_m'; % distance in # cells
            R_m = covariance_fcn(tau_m,len);
            [V,D] = eig(R_m);
            x_miu = zeros(size(tau));
            init_miu = x_miu;
            
            for p = 1:pf.Np
                %^^^^^^^^^^^^^^ initial noise ^^^^^^^^^^^^^^^^^^^^^^^
                sum_init_1 = 0; sum_init_2 = 0;
                for i = 1:size(V,1)
                    sum_init_1 = sum_init_1 + sqrt(D(i,i)) * randn * pf.init_stdev * V(:,i)';
                    sum_init_2 = sum_init_2 + sqrt(D(i,i)) * randn * pf.init_stdev * V(:,i)';
                end
                noise_init_1 = init_miu + sum_init_1; % a random field at time n
                noise_init_2 = init_miu + sum_init_2;
                x(:,:,p) = [U0_est(1,:) + noise_init_1; U0_est(2,:) + noise_init_2];
                x(x<0)=0;
                
            end
            
            % ******************* Uncorrelated **********************
        else
            for p = 1:pf.Np
                %^^^^^^^^^^^^^^ initial noise ^^^^^^^^^^^^^^^^^^^^^^^
                x(:,:,p) = U0_est + normrnd(0,pf.init_stdev,[size(U0_est,1), size(U0_est,2)]);
                x(x<0)=0;
                
            end
        end
% ================ Perform time propagation to obtain priori particles ================
        tic
        x_next = x; y_next = zeros(size(measure(U_est{1},pf))); N_eff = zeros(1,model_true.M);
        U_res{1} = U0_est;
    %  ****************************** PARTICLE FILTER STARTS ***********************************
        for n = 2:model_true.M
            %     fprintf('running time = %d\n',n);
            x(x<0)=0;
            d1l_temp = model_est.d1l(n);
            d2l_temp = model_est.d2l(n);
            d1r_temp = model_est.d1r(n);
            d2r_temp = model_est.d2r(n);
            
            for p = 1:pf.Np
                %^^^^^^^^^^^^^^ process noise ^^^^^^^^^^^^^^^^^^^^^^^^^^^
                % ///////////////////////////////////////////////////////
                if spatial_correlation
                    sum_x_1 = 0; sum_x_2 = 0; sum_y = 0;
                    for i = 1:size(V,1)
                        sum_x_1 = sum_x_1 + sqrt(D(i,i)) * randn * pf.model_stdev(1) * V(:,i)';
                        sum_x_2 = sum_x_2 + sqrt(D(i,i)) * randn * pf.model_stdev(2) * V(:,i)';
                    end
                    noise_x_1 = x_miu + sum_x_1;
                    noise_x_2 = x_miu + sum_x_2;
                end
                % ///////////////////////////////////////////////////////
                
                %^^^^^^^^^^^^^^ measurement noise ^^^^^^^^^^^^^^^^^^^^^^^
                noise_y = y_miu + normrnd(0,pf.meas_stdev,size(y_miu));
                
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
                flow_next = solver_flow(n-1,x(:,:,p),model_est);
                
                if spatial_correlation
                    x_next(:,:,p) = [flow_next(1,:) + noise_x_1; flow_next(2,:) + noise_x_2]; %correlated
                else
                    x_next(:,:,p) = flow_next + normrnd(0,pf.model_stdev(1),[size(flow_next,1),size(flow_next,2)]); %uncorrelated
                end
                x_next(x_next<0)=0;
                
                %------ measurement update ---------
                y_next(:,:,p) = measure(x_next(:,:,p),pf) + noise_y;
                y_next(y_next<0) = 0;
                y_true = [U_meas_true{n}(1,:)';U_meas_true{n}(2,:)']; % true measurement
                h = [y_next(1,:,p)'; y_next(2,:,p)']; % measurement equation
                m = size(h,1)*size(h,2);
                wt(p) = (2 * pi)^(-m/2) * (sqrt(sum(sum(abs(R).^2))))^(-1/2) * ...
                    exp(-1/2 * (y_true - h)'* R^(-1) * (y_true - h));
            end
            
            %------ normalize weight ----
            wt = wt./sum(wt); % make sure all the weights of each state i sum up to one
            N_eff(n) = (sum(wt.^2))^(-1); % effective particle size [HOW TO AVOID THE CURSE OF DIMENSIONALITY: SCALABILITY OF PARTICLE FILTERS WITH AND WITHOUT IMPORTANCE WEIGHTS?]
            %     [~,ind] = max(wt); % MAP estimate
            %     disp(N_eff(n))
            
            % ------ resampling -----
            for p = 1:pf.Np
                x_next(:,:,p) = x_next(:,:,(find(rand <= cumsum(wt),1)));
            end
            %   ******************************* PARTICLE FILTER ENDS **********************************
            
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
%                     drawnow
                    filename = sprintf('pf_test%d_%03d',test,n);
                    path = fullfile(directory,filename);
                    saveas(gca,path,'png')
                end
            end
            
            x = x_next;
        end
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
        
        save(sprintf('test%d_run%d',test,r));
        
    end % end of all runs
    log_beta_c1 = sum(log10(beta_c1))/Nr;
    log_beta_c2 = sum(log10(beta_c2))/Nr;
    true_beta_c1 = 10^log_beta_c1;
    true_beta_c2 = 10^log_beta_c2;
    
    %% save results to a table - influence of particles
    T_pi = [T_pi; Npc(pc) mean(runtime) true_beta_c1 true_beta_c2 mean(N_eff_mean) mean(MAE(5,:)) mean(MAE(6,:))];
    % clear MAE runtime N_eff_mean
    
end % end of all particle choices
colNames = {'particles','runtime','beta_1','beta_2','mean_N_eff','improvement_1','improvement_2'};
T_pi = array2table(T_pi,'VariableNames',colNames);

