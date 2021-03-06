%ASMSNMPCDISTCOLA Summary of this function goes here
% 
% [OUTPUTARGS] = ASMSNMPCDISTCOLA(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2019/01/15 18:03:46 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2019


global N scenario;
% number of mpc iteration
mpciterations = 150;
% number of prediction horizon
N             = 30;
% sampling time
T             = 1;  % [minute]
% initial controls (different initial conditions)
load Xinit29.mat;  
u0            = Xinit29(85:89);
u0            = repmat(u0,1,N);
% get initial measurement (states) at time T = 0.
tmeasure      = 0.0;
xmeasure      = Xinit29(1:84);

% number of realization and robust horizon
scr.numPar  = 2; % from +- 10% of price information
scr.nR      = 2; % robust horizon
scr.numScr  = 2^2;

% screnario tree NMPC dedicated 
[~, xmAs, umAs, obj, runtime] = sPfNmpcDistColA(@optProblem, @system, mpciterations, N, T, tmeasure, xmeasure, u0, scr);

%save MsAsNMPC21062019.mat;
%save MsAsNMPC04082019.mat;
%save MsAsNMPC03092019.mat;
keyboard;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the NMPC functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = system(t, x, u, T)
    
    global uc;
    uc = u;
    [~,x_out] = ode15s('cola_lv_cstr',[t t+T], x);
    lengthx   = size(x_out); 
    y         = x_out(lengthx(1),:)'; 
    
end

function [J,g,w0,w,lbg,ubg,lbw,ubw,params] = optProblem(x, u, N, x0_measure)   %add prediction horizon 
    import casadi.*
    
    % the model
    NT = 41;
    Uf = 0.3;           % Feeding rate F_0
    
    % invoke the model
    [~,state,xdot,inputs] = DistColACstr(Uf);
    f = Function('f',{state,inputs}, {xdot});
    
    % bound constraints
    VB_max = 4.008;
    
    % State bounds and initial guess
    x_min =  zeros(84,1);  % try without epsilon here, later put epsilon
    x_max =  ones(84,1);

    % Control bounds
    u_min = [0.1; 0.1; 0.1; 0.1; 0.1];
    u_max = [10; VB_max; 10; 1.0; 1.0];
    
    % compact bound variable for ease of function invocation 
    params.bound.x_min = x_min;
    params.bound.x_max = x_max;
    params.bound.u_min = u_min;
    params.bound.u_max = u_max;
    
    % Construct objective function
    load CstrDistXinit.mat;
    xf    = Xinit(1:84);
    u_opt = Xinit(85:89);
    
    % prices
    global mpciter;
    pf = 1; 
    pV = 0.02;
    pB = 2; 
    pD = 0;
    
    if mpciter == 50
        pf = 2;
        pV = 0.01;
        pB = 1;
        pD = 0;
    end
    
    if mpciter == 100
        pf = 4;
        pV = 0.05;
        pB = 3;
        pD = 0;
    end

    
    % compact price variable
    params.price.pf = pf;
    params.price.pV = pV;
    params.price.pB = pB;
    params.price.pD = pD;
    params.price.F_0= Uf;
    
    
    % dimensions
    global nx nu nk d tf ns;
    nx = 84;   % CSTR + Distillation Column A
    nu = 5;    % LT, VB, F, D, B
    nk = 1;
    tf = 1;      % in [minutes]
    h  = tf/nk;
    ns = 0;
    
    % compact model variable
    params.model.NT = NT;
    params.model.f  = f;
    params.model.xdot_val_rf_ss = xf;
    params.model.x  = x;
    params.model.u_opt = u_opt;
    % duplicate u0
    params.model.u  = repmat(u,1,nk); 


    % preparing collocation matrices
    [~,C,D,d] = collocationSetup();
    
    % compact collocation variable
    params.colloc.C = C;
    params.colloc.D = D;
    params.colloc.h = h;

    % start with an empty NLP
    w   = {};      % decision variables contain both control and state variables
    w0  = [];      % initial guess
    lbw = [];      % lower bound for decision variable
    ubw = [];      % upper bound
    J   = 0;       % objective function
    g   = {};      % nonlinear constraint
    lbg = [];      % lower bound for nonlinear constraint
    ubg = [];      % upper bound

    %delta_time = 60; % [minute] convert second to minute
    delta_time = 1;
    alpha = 1;
    beta  = 1;
    gamma = 1;
    
    % compact weight variable
    params.weight.delta_time = delta_time;
    params.weight.alpha      = alpha;
    params.weight.beta       = beta;
    params.weight.gamma      = gamma;
    
    % "Lift" initial conditions
    X0  = SX.sym('X0', nx);
    w   = {w{:}, X0};
    lbw = [lbw; x_min];
    ubw = [ubw; x_max];
    w0  = [w0; x(1,1:nx)'];
    g   = {g{:}, X0 - x0_measure};  % USE MEASUREMENT HERE !
    lbg = [lbg; zeros(nx,1)];
    ubg = [ubg; zeros(nx,1)];
    
    
    % Steady-state constraint
    global scenario;
    [ceq,xs] = buildSteadyStateConstraintsCstrDisColA(scenario);
    g        = {g{:}, ceq{:}};
    w        = {w{:}, xs{:}};
    w0       = [w0; x0_measure; u(:,1)];
    
    % bound constraints
    lb_u = [0.1; 0.1; 0.1; 0.1; 0.1];
    ub_u = [10; 4.008; 10; 1.0; 1.0];
    
    % State bounds and initial guess
    x_min     = zeros(84,1);
    x_max     = ones(84,1);
    xB_max    = 0.1;
    x_max(1)  = xB_max;
    x_min(84) = 0.3;
    x_max(84) = 0.7;
    lbx       = [x_min;lb_u];
    ubx       = [x_max;ub_u];
    NT        = 41;
    lbg       = [lbg;zeros(2*NT+2,1)];
    ubg       = [ubg;zeros(2*NT+2,1)];
    lbw       = [lbw; lbx];
    ubw       = [ubw; ubx];
    
    % formulate the NLP
    Xk = X0;

    
    load Qmax.mat;
    params.Qmax = Qmax;
    
    % HERE SHOULD BE LOOP N-TIMES ACCORDING TO THE NUMBER OF PREDICTION HORIZON
    count  = 2; % counter for state variable as initial guess
    ssoftc = 0;
    for i=1:N
        [J,g,w0,w,lbg,ubg,lbw,ubw,Xk,params,count,ssoftc] = iterateOnPredictionHorizon(Xk, w, w0, lbw, ubw, lbg, ubg, g, J, params, i, count,ssoftc);
    end

    
end


function [J,g,w0,w,lbg,ubg,lbw,ubw,Xk,params,count,ssoftc] = iterateOnPredictionHorizon(Xk, w, w0, lbw, ubw, lbg, ubg, g, J, params, iter, count, ssoftc)

   import casadi.*
   global N;
   % extract compact variables
   x_min = params.bound.x_min;
   x_max = params.bound.x_max;
   u_min = params.bound.u_min;
   u_max = params.bound.u_max;
   
   NT = params.model.NT;
   f  = params.model.f;
   xdot_val_rf_ss = params.model.xdot_val_rf_ss;
   x = params.model.x;
   u = params.model.u;
   u_opt = params.model.u_opt;
   
   pf = params.price.pf;
   pV = params.price.pV;
   pB = params.price.pB;
   pD = params.price.pD;
   F_0= params.price.F_0;
   
   
   C = params.colloc.C;
   D = params.colloc.D;
   h = params.colloc.h;
   
   delta_time = params.weight.delta_time;
   Qmax = params.Qmax;
   
   global nx nu nk d ns;

   for k=0:nk-1
        % New NLP variable for the control
        %Uk  = MX.sym(['U_' num2str(k)], nu);
        Uk     = SX.sym(['U_' num2str((iter-1)*nk+k)], nu);
        w      = {w{:}, Uk};
        lbw    = [lbw; u_min];
        ubw    = [ubw; u_max];
        indexU = (iter-1)*nk + (k+1);
        w0     = [w0;  u(:,indexU)];
        
        Jcontrol   = (Qmax(nx+1:nx+nu,1).*(Uk - u_opt))' * (Uk - u_opt);
      
        % State at collocation points
        Xkj   = {};
        SumX1 = 0;
        for j=1:d
            Xkj{j} = SX.sym(['X_' num2str((iter-1)*nk+k) '_' num2str(j)], nx);
            w      = {w{:}, Xkj{j}};
            lbw    = [lbw; x_min];
            ubw    = [ubw; x_max];   
            w0     = [w0; x(iter+1,:)'];  
            count  = count + 1;
        end

        % Loop over collocation points
        Xk_end = D(1)*Xk; 

        for j=1:d
           % Expression for the state derivative at the collocation point
           xp = C(1,j+1)*Xk;
           for r=1:d
               xp = xp + C(r+1,j+1)*Xkj{r};
           end

           % Append collocation equations
           fj  = f(Xkj{j},Uk);
           g   = {g{:}, h*fj - xp};
           lbg = [lbg; zeros(nx,1)];
           ubg = [ubg; zeros(nx,1)];

           % Add contribution to the end state
           Xk_end = Xk_end + D(j+1)*Xkj{j};
           
        end    

        % New NLP variable for state at end of interval
        Xk  = SX.sym(['X_' num2str((iter-1)*nk+k)], nx);
        w   = {w{:}, Xk};
        lbw = [lbw; x_min];
        ubw = [ubw; x_max];
        w0  = [w0; x(iter+1,:)'];
        count  = count + 1;

        % Add equality constraint
        g   = {g{:}, (Xk_end-Xk)};
        lbg = [lbg; zeros(nx,1)];
        ubg = [ubg; zeros(nx,1)];
               
        
        Jstate =(Qmax(1:nx,1).*(Xk - xdot_val_rf_ss))' * (Xk - xdot_val_rf_ss) * delta_time;
        
        global scenario;
        switch scenario
            case 1
                Jecon  = ((pf+0.5*pf)*F_0 + (pV+0.5*pV)*Uk(2) - (pB+0.5*pB)*Uk(5) - pD*Uk(4)) * delta_time;
            case 2
                Jecon  = ((pf-0.5*pf)*F_0 + (pV-0.5*pV)*Uk(2) - (pB-0.5*pB)*Uk(5) - pD*Uk(4)) * delta_time;
            case 3
                Jecon  = ((pf+0.25*pf)*F_0 + (pV+0.25*pV)*Uk(2) - (pB+0.25*pB)*Uk(5) - pD*Uk(4)) * delta_time;
            case 4
                Jecon  = ((pf-0.25*pf)*F_0 + (pV-0.25*pV)*Uk(2) - (pB-0.25*pB)*Uk(5) - pD*Uk(4)) * delta_time;
        end
        
        alpha  = 1;
        beta   = 1;
        gamma  = 1;
        
        J = J + alpha*Jcontrol + gamma*Jstate + beta*Jecon;

    end
end