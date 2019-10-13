function U_n = solver_flow(n,U,model_param)
% propagate the flow (state vector) to the next timestep given initial
% conditions and boundary conditions

% the boundary condition
if contains(model_param.model_name,'LWR')==false % if not LWR model
    lbc = [model_param.d1l(n+1);model_param.d2l(n+1)];
    rbc = [model_param.d1r(n+1);model_param.d2r(n+1)];
    % insert the boundary condition
    Up1=[U(:,2:end),rbc];    % downstream initial states
    Um1=[lbc,U(:,1:end-1)];  % upstream initial states
end

% original code by Shimao Fan can be found https://github.com/shimaof/heterogeneous-traffic-model
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% phase transition happens (if creeping model)
%////////////////////////////

%============================================
% updating with godunov methods/heterogeneous extension of Cell
% Transmission Model

switch model_param.model_name
    case {'Creeping model','True model'}
        rl = Um1(1,:) + Um1(2,:);   % total density on left
        r = U(1,:) + U(2,:);        % total density in middle
        rr = Up1(1,:) + Up1(2,:);   % total density on right
        
        case1 = rl<=model_param.rm2;        % left state in the domain D1
        case2 = r<=model_param.rm2;         % middle state in the domain D1
        case3 = rr<=model_param.rm2;        % right state in the domain D1
        
        id_l = find(case1-case2); % find all phase changes ul and u
        id_r = find(case2-case3); % find all phase changes u and ur
        
        idl = case1(id_l)-case2(id_l);% -1 means left large density
        idr = case2(id_r)-case3(id_r);% +1 mean left state small
        
        for i = 1:length(idr)-1
            if idr(i)<0
                U(1,id_r(i)) = model_param.rm2-U(2,id_r(i));   % insert a intermediate state
            else
                Up1(1,id_r(i)) = model_param.rm2-Up1(2,id_r(i));   % insert a intermediate state
            end
        end
        
        U_n = U+model_param.lambda*(Num_Flux(Um1,U,model_param)- ...
            Num_Flux(U,Up1,model_param));
    case 'n-populations model'
        U_n = U+model_param.lambda*(Num_Flux(Um1,U,model_param)- ...
            Num_Flux(U,Up1,model_param));
    case 'Porous model'
        U_n = U+model_param.lambda*(HLL_porous(Um1,U,model_param.vm,model_param.s1,model_param.s2,model_param.lda,model_param.c1,model_param.c2,model_param.alpha1,model_param.alpha2)...
            -HLL_porous(U,Up1,model_param.vm,model_param.s1,model_param.s2,model_param.lda,model_param.c1,model_param.c2,model_param.alpha1,model_param.alpha2));
    case 'LWR'
        U_n = U;
        U1 = [model_param.d1l(n), U(1,:), model_param.d1r(n)];
        U2 = [model_param.d2l(n), U(2,:), model_param.d2r(n)];
        for i = 2:length(U1)-1
            U1_n(i-1) = U1(i) + model_param.dt/model_param.dx * (min(sending_greenshield(U1(i-1),model_param.vm1,model_param.rm1),...
            receiving_greenshield(U1(i),model_param.vm1,model_param.rm1)) - min(sending_greenshield(U1(i),model_param.vm1,model_param.rm1), ...
            receiving_greenshield(U1(i+1),model_param.vm1,model_param.rm1)));
            
            U2_n(i-1) = U2(i) + model_param.dt/model_param.dx * (min(sending_greenshield(U2(i-1),model_param.vm2,model_param.rm2),...
            receiving_greenshield(U2(i),model_param.vm2,model_param.rm2)) - min(sending_greenshield(U2(i),model_param.vm2,model_param.rm2), ...
            receiving_greenshield(U2(i+1),model_param.vm2,model_param.rm2)));
        end
        U_n = [U1_n; U2_n];
end
%============================================

% end


%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% subfunctionals
%///////////////////////////////////////////////////////
function y = new_flux(Ul,Ur,model_param,num_param)
rl = Ul(1,:) + Ul(2,:); % total density on left
rr = Ur(1,:) + Ur(2,:); % total density on right
%     n = length()
inst1 = (rl==0 | rr==0); % case either initial data are zero
inst2 = 1-inst1;   % other cases

% update the numerical flux
flux_hll = HLL(Ul,Ur,model_param,num_param);      % case inst2
flux_sr = Num_Flux(Ul,Ur,model_param);  % sending receiving
indx1 = find(isnan(flux_hll(1,:)));
indx2 = find(isnan(flux_hll(2,:)));
flux_hll(1,indx1) = flux_sr(1,indx1);
flux_hll(2,indx2) = flux_sr(2,indx2);
%     y = flux_hll;
y = (ones(2,1)*inst1).*flux_sr + (ones(2,1)*inst2).*flux_hll;
function F = HLL(Ul,Ur,model_param,num_param) % the HLL Riemann solver

% where v contains boundary condition of velocity
% eigenvalues:
eigns1 = lambda(Ul,model_param,num_param);  % left bound
eigns2 = lambda(Ur,model_param,num_param);  % right bound
% another way to define c1 and c2
c1 = min(eigns1(1,:),eigns2(1,:));
c2 = max(eigns1(2,:),eigns2(2,:));
case1 = c1>=0;
case2 = c1<0 & c2>=0;
case3 = c2<0;
% left flux
Fl = flux(Ul,model_param);  % plug into left boundary velocity information
% right flux
Fr = flux(Ur,model_param);
%-----------------------------------------------
% intermediate flux
F_hll(1,:) = (c2.*Fl(1,:)-c1.*Fr(1,:)+(Ur(1,:)-Ul(1,:)).*c1.*c2)./...
    (c2-c1);
F_hll(2,:) = (c2.*Fl(2,:)-c1.*Fr(2,:)+(Ur(2,:)-Ul(2,:)).*c1.*c2)./...
    (c2-c1);
F(1,:) = case1.*Fl(1,:)+case2.*F_hll(1,:)+case3.*Fr(1,:);
F(2,:) = case1.*Fl(2,:)+case2.*F_hll(2,:)+case3.*Fr(2,:);

%--------------------------------------------------------------------------
function y = lambda(U,model_param)  % U = (density, w), calculate eigenvalues

r = U(1,:)+U(2,:);    % sum
u1 = vel(r,model_param.vm1,model_param.rm1);
u2 = max(vel(r,model_param.vm2,model_param.rm2),0);
dv1 = DiffVel(r,model_param.vm1,model_param.rm1);
dv2 = DiffVel(r,model_param.vm2,model_param.rm2);
indx = find(r>model_param.rm2);
dv2(indx) = 0;
% characteristic speed
c1 = U(1,:).*dv1;
c2 = U(2,:).*dv2;
q1 = c1+u1; q2 = u2+c2;
delta = sqrt((q1-q2).^2+4*c1.*c2);
y(1,:) = .5*(q1+q2-delta);   % slower chareateritic field
y(2,:) = .5*(q1+q2+delta);   % faster characteristic field

%--------------------------------------------------------------
function  F = Num_Flux(Ul,Ur,model_param) % numerical flux
% find sending function/first class
s_1 = sending_1(Ul,model_param);
% receiving of flux/first class
r_1 = receiving_1(Ul,Ur,model_param);

% find sending function
s_2 = sending_2(Ul,model_param);
% receiving of flux
r_2 = receiving_2(Ul,Ur,model_param);
% conservation of mass
F(1,:) = min(s_1,r_1);   % take the minimum, flow of first class
% conservation of momentum
F(2,:) = min(s_2,r_2);   % take the minimum, flow of the second class
%-------------------------------------------------
% define the receiving function for 1st vehicle class
function y = receiving_1(Ul,Ur,model_param)

rhoc = (model_param.rm1-Ur(2,:))/2;  % critical density
case1 = Ur(1,:)<=rhoc;  % maximum possible
case2 = Ur(1,:)>rhoc;   % maximum
u0(1,:) = rhoc;
u0(2,:) = Ur(2,:);
q_max = traffic_flux_1(u0,model_param); % maximum flow rate
%       q_max = max(traffic_flux_1(Ur,model_param.vm1,model_param.rm1)); % maximum flow rate
q = traffic_flux_1(Ur,model_param);
y = case2.*q + case1.*q_max;
%--------------------------------------------------
% define the receiving function for 2nd vehicle class
function y = receiving_2(Ul,Ur,model_param)
% when r<model_param.rm2, non-creeping phase
r = Ur(1,:) + Ur(2,:); % total density

rhoc = max((model_param.rm2-Ur(1,:))/2,0);  % critical density
case1 = Ur(2,:)<=rhoc;  % maximum possible
case2 = Ur(2,:)>rhoc;   % maximum
u0(2,:) = rhoc;
u0(1,:) = Ur(1,:);
q_max = traffic_flux_2(u0,model_param); % maximum flow rate
q = traffic_flux_2(Ur,model_param);
y = case2.*q + case1.*q_max;
%----------------------------------------------
% define sending function for 1st vehicle class
function y = sending_1(U,model_param)% sending of creeping vehicles/for two phases

rhoc = (model_param.rm1-U(2,:))/2;  % a vector of rhoc/3-params fd
case1 = U(1,:)<=rhoc;  % maximum possible
case2 = U(1,:)>rhoc;   % maximum
u0(1,:) = rhoc;
u0(2,:) = U(2,:);
q_max = traffic_flux_1(u0,model_param); % maximum flow rate
%       q_max = max(traffic_flux_1(U,model_param.vm1,model_param.rm1));
q = traffic_flux_1(U,model_param);
y = case1.*q + case2.*q_max;
%------------------------------------------------
% define sending function for 2nd vehicle class
function y = sending_2(U,model_param)% U = [\rho, w] two state variables
rhoc = max((model_param.rm2-U(1,:))/2,0);  % critical density
case1 = U(2,:)<=rhoc;    % maximum possible
case2 = U(2,:)>rhoc;     % maximum
u0(2,:) = rhoc;
u0(1,:) = U(1,:);
q_max = traffic_flux_2(u0,model_param); % maximum flow rate
q = traffic_flux_2(U,model_param);
y = case1.*q + case2.*q_max;
%-----------------------------------------------
function y = flux(U,model_param)  % define flux function
r = U(1,:)+U(2,:);
u1 = vel(r,model_param.vm1,model_param.rm1);
u2 = vel(r,model_param.vm2,model_param.rm2);
y(1,:) = U(1,:).*u1;
y(2,:) = U(2,:).*u2;

%-------------------------------------------------------
function y = traffic_flux_1(U,model_param)  % traffic flux of 1st class
rho = U(1,:)+U(2,:);
v = vel(rho,model_param.vm1,model_param.rm1);
y = U(1,:).*v;
%-------------------------------------------------------
function y = traffic_flux_2(U,model_param) % traffic flux of 2nd class
rho = U(1,:)+U(2,:);
v = vel(rho,model_param.vm2,model_param.rm2);
y = U(2,:).*v;
y = max(y,0*rho);
%-------------------------------------------------------
function y = DiffVel(r,vm,rm)% define derivative of velocity
y = 0*r-vm./rm;
%-------------------------------------------------------
function y = vel(rho,vm,rhom)  % define velocity functions
y = vm*(1-(rho/rhom));    % which is the quadratic form
%-------------------------------------------------------










function F = HLL_porous(Ul,Ur,vm,s1,s2,lda,k1,k2,alpha1,alpha2) % the HLL Riemann solver
% where v contains boundary condition of velocity

% eigenvalues:
eigns1 = lambda_porous(Ul,vm,s1,s2,lda,k1,k2,alpha1,alpha2);  % left bound
eigns2 = lambda_porous(Ur,vm,s1,s2,lda,k1,k2,alpha1,alpha2);  % right bound
% another way to define c1 and c2
c1 = min(eigns1(1,:),eigns2(1,:));
c2 = max(eigns1(2,:),eigns2(2,:));

case1 = c1>=0;
case2 = c1<0 & c2>=0;
case3 = c2<0;
case4 = (c1==c2);

Fl = max(flux_porous(Ul,vm,s1,s2,lda,k1,k2,alpha1,alpha2),0);  % plug into left boundary velocity information
Fr = max(flux_porous(Ur,vm,s1,s2,lda,k1,k2,alpha1,alpha2),0);
%-----------------------------------------------
% intermediate flux
F_hll(1,:) = (1-case4).*(c2.*Fl(1,:)-c1.*Fr(1,:)+(Ur(1,:)-Ul(1,:)).*c1.*c2)./...
    (c2-c1) + case4*0;

F_hll(2,:) = (1-case4).*(c2.*Fl(2,:)-c1.*Fr(2,:)+(Ur(2,:)-Ul(2,:)).*c1.*c2)./...
    (c2-c1)+case4*0;
F(1,:) = case1.*Fl(1,:)+case2.*F_hll(1,:)+case3.*Fr(1,:);
F(2,:) = case1.*Fl(2,:)+case2.*F_hll(2,:)+case3.*Fr(2,:);


%--------------------------------------------------------------------------
function y = lambda_porous(U,vm,s1,s2,lambda,c1,c2,alpha1,alpha2)  % U = (density, w), calculate eigenvalues
r = U(1,:)+U(2,:);    % sum
%        u1 = max(vel(r,vm,s1,lambda,c,alpha1,alpha2),0);
%        u2 = max(vel(r,vm,s2,lambda,c,alpha1,alpha2),0);
u1 = vel_porous(r,vm,s1,lambda,c1,alpha1,alpha2);
u2 = vel_porous(r,vm,s2,lambda,c2,alpha1,alpha2);
dv1 = DiffVel_porous(r,vm,s1,lambda,c1,alpha1,alpha2);
dv2 = DiffVel_porous(r,vm,s2,lambda,c2,alpha1,alpha2);
%        indx = find(r>=rm2);
%        dv2 = DiffVel(r,vm2,rm2);
%        dv2(indx) = 0;
% characteristic speed
c1 = U(1,:).*dv1;
c2 = U(2,:).*dv2;
q1 = c1+u1; q2 = u2+c2;
delta = sqrt((q1-q2).^2+4*c1.*c2);

y(1,:) = .5*(q1+q2-delta);   % slower chareateritic field
y(2,:) = .5*(q1+q2+delta);   % faster characteristic field
%-----------------------------------------------
function y = flux_porous(U,vm,s1,s2,lambda,c1,c2,alpha1,alpha2)  % define flux function
r = U(1,:)+U(2,:);
u1 = vel_porous(r,vm,s1,lambda,c1,alpha1,alpha2);
u2 = vel_porous(r,vm,s2,lambda,c2,alpha1,alpha2);
y(1,:) = U(1,:).*u1;
y(2,:) = U(2,:).*u2;
%-----------------------------------------------
function y = DiffVel_porous(r,vm,sj,lambda,c,alpha1,alpha2)% define derivative of velocity
eps = 1.0e-03;
y = (vel_porous(r+eps,vm,sj,lambda,c,alpha1,alpha2)-vel_porous(r,vm,sj,lambda,c,alpha1,alpha2))/eps;

%===============================================
function y = vel_porous(r,vm,sj,lambda,c,alpha1,alpha2)
vf = v_porous(r,vm,sj,lambda,c,alpha1);
vs = v_porous(r,vm,sj,lambda,c,alpha2);

y1 = (1-sj./lambda.*exp(c.*r)).*vf + sj./lambda.*exp(c.*r).*vs;
y = max(y1,0);

function y = v_porous(r,vm,sj,lambda,c,alpha)

y = vm.*(1-sj./lambda.*exp(c.*r)).^alpha;




function [send_flow] = sending_greenshield(rho, vm,rm)
   rho_c = rm/2;
   
   if rho <= rho_c
       send_flow = greenshield(rho,vm,rm);
   else
       send_flow = vm * rm/4;
   end
%-----------------------------------------------
function [receive_flow] = receiving_greenshield(rho, vm,rm)
rho_c = rm/2;
   if 0 <= rho && rho <= rho_c
       receive_flow = vm * rm/4;
   elseif rho_c < rho && rho <= rm
       receive_flow = greenshield(rho,vm,rm);
   else
       receive_flow = 0;
   end
   
function q = greenshield(rho,vm,rm)
if rho<0
    q = 0;
else
    q = vm*(1-rho/rm)*rho;
end