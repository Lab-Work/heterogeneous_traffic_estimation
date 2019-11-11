function U = initialize(model_param,test)
% initialize model parameters, discretization, initial and boundary
% conditions of model_param and test = 1,2,3 or 4.
% all the quantities are unit-less, i.e., normalized s.t. scalable
name = model_param.model_name;

switch name
    case 'True model' % true model is creeping model
        model_param.rm1 = 1.8; model_param.rm2 =1.0;
        model_param.vm1 = 1.8; model_param.vm2 = model_param.vm1;
        
    case 'Creeping model'
        model_param.rm1 = 1.7; model_param.rm2 = 0.9;
        model_param.vm1 = 1.9; model_param.vm2 = model_param.vm1;

    case 'Porous model'
        model_param.rm1 = 1.8; model_param.rm2 = 1.0;
        model_param.vm = 1.8;
        model_param.s1 = 1.0; model_param.s2 = 1.4; % critical pore size

        model_param.alpha1 = 1; model_param.alpha2 = 1;
        model_param.lda = 4;
        model_param.c1 = log(model_param.lda/model_param.s1)/model_param.rm1;
        model_param.c2 = log(model_param.lda/model_param.s2)/model_param.rm2;
end
model_param.tfinal = 3; % total simulation time
model_param.len = 2;
model_param.x = linspace(0,model_param.len,40);
model_param.dx = model_param.x(2)-model_param.x(1);
model_param.dt = (model_param.x(2)-model_param.x(1))/2.16;
model_param.lambda = model_param.dt/model_param.dx;
model_param.t = 0:model_param.dt:model_param.tfinal;
model_param.M = length(model_param.t);
model_param.N = length(model_param.x);

%==========================================
switch test
    case 1   % overtaking in freeflow
        thres1 = 2;
        thres2 = 0.4;
        thres3 = 0.8;
        freq = 0.07;
        switch name
            case 'True model'
                % Initial densities
                U(1,1:thres1) = eps;
                U(1,thres1+1:ceil(thres2*model_param.N/model_param.len)) = .5;
                U(1,ceil(thres2*model_param.N/model_param.len)+1:model_param.N) = eps;
                U(2,1:ceil(thres2*model_param.N/model_param.len)) = eps;
                U(2,ceil(thres2*model_param.N/model_param.len)+1:ceil(thres3*model_param.N/model_param.len)) = .6;
                U(2,ceil(thres3*model_param.N/model_param.len)+1:model_param.N) = eps;
                % Boundary densities
                model_param.d1l = (square(freq*[1:model_param.M])*0.04+0.1)'; % upstream
                model_param.d2l = (square(freq*[1:model_param.M])*0.04+0.1)'; 
                model_param.d1r = 0+0*model_param.t; % downstream
                model_param.d2r = 0+0*model_param.t;
                
            case 'Creeping model'
                U(1,1:thres1) = eps;
                U(1,thres1+1:ceil(thres2*model_param.N/model_param.len)) = .7;
                U(1,ceil(thres2*model_param.N/model_param.len)+1:model_param.N) = eps;
                U(2,1:ceil(thres2*model_param.N/model_param.len)) = eps;
                U(2,ceil(thres2*model_param.N/model_param.len)+1:ceil(thres3*model_param.N/model_param.len)) = .7;
                U(2,ceil(thres3*model_param.N/model_param.len)+1:model_param.N) = eps;

                model_param.d1l = (square(freq*[1:model_param.M])*0.04+0.04)';
                model_param.d2l = (square(freq*[1:model_param.M])*0.04+0.04)';
                model_param.d1r = 0.1+0*model_param.t;
                model_param.d2r = 0.1+0*model_param.t;
                
            case 'Porous model'
                U(1,1:thres1) = eps;
                U(1,thres1+1:ceil(thres2*model_param.N/model_param.len)) = .5;
                U(1,ceil(thres2*model_param.N/model_param.len)+1:model_param.N) = eps;
                U(2,1:ceil(thres2*model_param.N/model_param.len)) = eps;
                U(2,ceil(thres2*model_param.N/model_param.len)+1:ceil(thres3*model_param.N/model_param.len)) = .7;
                U(2,ceil(thres3*model_param.N/model_param.len)+1:model_param.N) = eps;
                
                model_param.d1l = (square(0.08*[1:model_param.M])*0.08+0.08)';
                model_param.d2l = (square(0.08*[1:model_param.M])*0.08+0.08)';
                model_param.d1r = U(1,end)+0*model_param.t;
                model_param.d2r = U(2,end)+0*model_param.t;
                
        end
       
    case 2  % congested traffic
        thres1 = 15;
        thres2 = 22;
        freq = 0.07;
        switch name  
            case 'True model'
                U(1,1:thres1) = 0.0;
                U(1,thres1+1:thres2) = 0.55;
                U(1,thres2+1:model_param.N) = 0.0;
                U(2,:) = 0.9;

                model_param.d1l = U(1,1)+0*model_param.t;
                model_param.d2l = U(2,1)+0*model_param.t;
                model_param.d1r = U(1,end)+0.1+0*model_param.t;
                model_param.d2r = U(2,end)+0*model_param.t;
                
            case 'Creeping model'
                U(1,1:thres1) = 0.1;
                U(1,thres1+1:thres2) = 0.45;
                U(1,thres2+1:model_param.N) = 0.1;
                U(2,:) = 0.7;

%                 model_param.d1l = U(1,1)+(square(freq*[1:model_param.M])*0.03+0.03)';
                model_param.d1l = U(1,1)+0*model_param.t;
                model_param.d2l = U(2,1)+0*model_param.t;
                model_param.d1r = U(1,end)+0*model_param.t;
                model_param.d2r = U(2,end)+0*model_param.t;
        end
  
     case 3  % queue clearance: congested to freeflow
        switch name
            case 'True model'
                U(1,1:ceil(model_param.N/2)+5) = 1.4*ones(1,ceil(model_param.N/2)+5);
                U(1,ceil(model_param.N/2)+5+1:model_param.N) = eps*ones(1,model_param.N-ceil(model_param.N/2)-5);
                U(2,1:ceil(model_param.N/2)+5) = 0.6*ones(1,ceil(model_param.N/2)+5);
                U(2,ceil(model_param.N/2)+5+1:model_param.N) = eps*ones(1,model_param.N-ceil(model_param.N/2)-5);

                model_param.d1l = U(1,1)+0*model_param.t;
                model_param.d2l = U(2,1)+0*model_param.t;
                model_param.d1r = 0.2+0*model_param.t;
                model_param.d2r = 0.2+0*model_param.t;
                
            case 'Creeping model'
                U(1,1:ceil(model_param.N/2)+5) = 1.2*ones(1,ceil(model_param.N/2)+5);
                U(1,ceil(model_param.N/2)+5+1:model_param.N) = eps*ones(1,model_param.N-ceil(model_param.N/2)-5);
                U(2,1:ceil(model_param.N/2)+5) = 0.4*ones(1,ceil(model_param.N/2)+5);
                U(2,ceil(model_param.N/2)+5+1:model_param.N) = eps*ones(1,model_param.N-ceil(model_param.N/2)-5);
                
                model_param.d1l = U(1,1)+0*model_param.t;
                model_param.d2l = U(2,1)+0*model_param.t;
                model_param.d1r = 0.1+0*model_param.t;
                model_param.d2r = 0.1+0*model_param.t;
        end

    case 4 % creeping traffic
        thres1 = 2;
        thres2 = 0.7;
        freq = 0.07;
        switch name
            case 'True model'
                U(1,1:thres1) = eps;
                U(1,thres1+1:ceil(thres2*model_param.N/model_param.len)) = .4;
                U(1,ceil(thres2*model_param.N/model_param.len)+1:model_param.N) = eps;
                U(2,1:ceil(thres2*model_param.N/model_param.len)) = eps;
                U(2,ceil(thres2*model_param.N/model_param.len)+1:model_param.N) = .8;
                
                model_param.d1l = (square(freq*[1:model_param.M])*0.04+0.04+0.04)';
                model_param.d2l = (square(freq*[1:model_param.M])*0.04+0.04+0.04)';
                model_param.d1r = 0.0+0*model_param.t;
                model_param.d2r = 1+0*model_param.t;
                
            case 'Creeping model'
                U(1,1:thres1) = eps;
                U(1,thres1+1:ceil(thres2*model_param.N/model_param.len)) = .3;
                U(1,ceil(thres2*model_param.N/model_param.len)+1:model_param.N) = eps;
                U(2,1:ceil(thres2*model_param.N/model_param.len)) = eps;
                U(2,ceil(thres2*model_param.N/model_param.len)+1:model_param.N) = .7;
                
                model_param.d1l = (square(freq*[1:model_param.M])*0.04+0.04)';
                model_param.d2l = (square(freq*[1:model_param.M])*0.04+0.04)';
                model_param.d1r = 0.1+0*model_param.t;
                model_param.d2r = 0.8+0*model_param.t;
                
            case 'Porous model'
                U(1,1:thres1) = eps;
                U(1,thres1+1:ceil(thres2*model_param.N/model_param.len)) = .5;
                U(1,ceil(thres2*model_param.N/model_param.len)+1:model_param.N) = eps;
                U(2,1:ceil(thres2*model_param.N/model_param.len)) = eps;
                U(2,ceil(thres2*model_param.N/model_param.len)+1:model_param.N) = .7;
                
                model_param.d1l = (square(freq*[1:model_param.M])*0.08+0.08)';
                model_param.d2l = (square(freq*[1:model_param.M])*0.08+0.08)';
                model_param.d1r = 0.0+0*model_param.t;
                model_param.d2r = 0.999+0*model_param.t;
        end
end
end