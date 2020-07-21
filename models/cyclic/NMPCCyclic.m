%function NMPCCyclic
%NMPCCYCLIC Summary of this function goes here
% 
% [OUTPUTARGS] = NMPCCYCLIC(INPUTARGS) Explain usage here
% 
% An economic NMPC controller for a nondissipative process.
%
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2018/01/17 21:35:17 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2018

import casadi.* 
format long;

global N delta le nx;
nx = 3;
le = ones(3,1);
% number of mpc iteration
%mpciterations = 40;
mpciterations = 100;
% number of prediction horizon
%N             = 50;
N             = 200;
%N             = 500;
%N             = 1000;
% sampling time
T             = 1;  % [minute]
% initial control and state
u0            = 0.3;
u0            = repmat(u0,1,N);
% get initial measurement (states) at time T = 0.
tmeasure      = 0.0;
%xmeasure      = [0.1;0.1;0.1];
xmeasure      = [1.0;1e-4;0.1];

% economic NMPC
[~, xmeasureAll, uAll, xuSS, optRes, params, runtime] = iNmpcCyclic(@optProblem, @system, mpciterations, N, T, tmeasure, xmeasure, u0);

keyboard;
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the NMPC functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = system(t, x, u, T)
    
    global uc;
    uc = u;
    [~,x_out] = ode15s('cyclic_ode',[t t+T], x);
    lengthx   = size(x_out); 
    y         = x_out(lengthx(1),:)'; 
    
end

function [J,g,w0,w,lbg,ubg,lbw,ubw,params] = optProblem(x, u, N, x0_measure, xuSS)   %add prediction horizon 
    import casadi.*
    
    % invoke the model
    [~,state,xdot,inputs] = cyclicFunc();
    f = Function('f',{state,inputs}, {xdot});
    
    % State bounds and initial guess
    x_min =  zeros(3,1);  
    x_max =  ones(3,1);
    
    % Control bounds
    u_min = 0.049;
    u_max = 0.449;
    
    % compact bound variable for ease of function invocation 
    params.bound.x_min = x_min;
    params.bound.x_max = x_max;
    params.bound.u_min = u_min;
    params.bound.u_max = u_max;
    
    % dimensions
    global nx nu nk d tf ns;
    nx = 3;   
    nu = 1;    
    nk = 1;
    tf = 0.1;   % in [minutes]
    %tf = 3;    % 3 minutes
    h  = tf/nk;
    ns = 0;
    
    % compact model variable
    params.model.f  = f;
    %params.model.xdot_val_rf_ss = xf;
    params.model.x  = x;
    %params.model.u_opt = u_opt;
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
    
    % incorporate steady-state optimization
    u1      = MX.sym('u1');     % u  - heat flux through cooling jacket
    x1      = MX.sym('x1');     % x1 - concentration of R
    x2      = MX.sym('x2');     % x2 - concentration of P1
    x3      = MX.sym('x3');     % x3 - concentration of P2 
    %epsilon = MX.sym('epsilon');
    xs     = [x1;x2;x3;u1]; 
    global delta;
    x1dot  = 1 - 1e4*x1^2*exp(-1/x3) - 400*x1*exp(-delta/x3) - x1;
    x2dot  = 1e4*x1^2*exp(-1/x3) - x2;
    x3dot  = u1 - x3;
    ceq    = [x1dot;x2dot;x3dot];
    g      = {g{:}, ceq};
    lbg    = [lbg; 0; 0; 0];
    ubg    = [ubg; 0; 0; 0];
    w      = {w{:}, xs};
    w0     = [w0; x0_measure; u(1)];
    lbx    = [0;0;0;0.049];
    ubx    = [1;1;1;0.449];
    lbw    = [lbw; lbx];
    ubw    = [ubw; ubx];
        
    % "Lift" initial conditions
    X0  = MX.sym('X0', nx);
    w   = {w{:}, X0};
    lbw = [lbw; x_min];
    ubw = [ubw; x_max];
    w0  = [w0; x(1,1:nx)'];
    g   = {g{:}, X0 - x0_measure};  % USE MEASUREMENT HERE !
    lbg = [lbg; zeros(nx,1)];
    ubg = [ubg; zeros(nx,1)];

    % formulate the NLP
    Xk = X0;
    
    % HERE SHOULD BE LOOP N-TIMES ACCORDING TO THE NUMBER OF PREDICTION HORIZON
    count  = 2; % counter for state variable as initial guess
    ssoftc = 0;
    for i=1:N
        [J,g,w0,w,lbg,ubg,lbw,ubw,Xk,params,count,ssoftc] = iterateOnPredictionHorizon(Xk, w, w0, lbw, ubw, lbg, ubg, g, J, params, i, count, ssoftc, xs);
    end

    
end


function [J,g,w0,w,lbg,ubg,lbw,ubw,Xk,params,count,ssoftc] = iterateOnPredictionHorizon(Xk, w, w0, lbw, ubw, lbg, ubg, g, J, params, iter, count, ssoftc, xs)

   import casadi.*
   global N;
   % extract compact variables
   x_min = params.bound.x_min;
   x_max = params.bound.x_max;
   u_min = params.bound.u_min;
   u_max = params.bound.u_max;
   
   f  = params.model.f;
   %xdot_val_rf_ss = params.model.xdot_val_rf_ss;
   x = params.model.x;
   u = params.model.u;
   %u_opt = params.model.u_opt;
  
   
   C = params.colloc.C;
   D = params.colloc.D;
   h = params.colloc.h;
   
   delta_time = params.weight.delta_time;
   %Qmax = params.Qmax;
   
   global nx nu nk d xf uf check Qc le;

   for k=0:nk-1
        % New NLP variable for the control
        %Uk  = MX.sym(['U_' num2str(k)], nu);
        Uk     = MX.sym(['U_' num2str((iter-1)*nk+k)], nu);
        w      = {w{:}, Uk};
        lbw    = [lbw; u_min];
        ubw    = [ubw; u_max];
        indexU = (iter-1)*nk + (k+1);
        w0     = [w0;  u(:,indexU)];
        
        %Jcontrol   = (Qmax(nx+1:nx+nu,1).*(Uk - u_opt))' * (Uk - u_opt);
        %Jcontrol   = 2*(Uk - 0.14909) * (Uk - 0.14909);
%         if check == 0
%             Jcontrol   = 2*(Uk - xs(4)) * (Uk - xs(4));
%         else
%             Jcontrol   = 2*(Uk - uf) * (Uk - uf);
%         end
%         if check == 0
%             Jcontrol   = 4*(Uk - xs(4)) * (Uk - xs(4));
%         else
%             Jcontrol   = 4*(Uk - uf) * (Uk - uf);
%         end
%         if check == 0
%             Jcontrol   = 10*(Uk - xs(4)) * (Uk - xs(4));
%         else
%             Jcontrol   = 10*(Uk - uf) * (Uk - uf);
%         end
%         if check == 0
%             Jcontrol   = 10*(Uk - xs(4)) * (Uk - xs(4));
%             %Jcontrol   = 0;
%         else
%             %Jcontrol   = Qc(4,1)*(Uk - uf) * (Uk - uf);
%             %Jcontrol   = 10*(Uk - xs(4)) * (Uk - xs(4));
%             Jcontrol   = 10*(Uk - uf) * (Uk - uf);
%         end

        Jcontrol   = 10*(Uk - xs(4)) * (Uk - xs(4));
        %Jcontrol   = 5*(Uk - xs(4)) * (Uk - xs(4));
      
        % State at collocation points
        Xkj   = {};
        SumX1 = 0;
        for j=1:d
            Xkj{j} = MX.sym(['X_' num2str((iter-1)*nk+k) '_' num2str(j)], nx);
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
        Xk  = MX.sym(['X_' num2str((iter-1)*nk+k)], nx);
        w   = {w{:}, Xk};
        lbw = [lbw; x_min];
        ubw = [ubw; x_max];
        w0  = [w0; x(iter+1,:)'];
        count  = count + 1;

        % Add equality constraint
        g   = {g{:}, (Xk_end-Xk)};
        lbg = [lbg; zeros(nx,1)];
        ubg = [ubg; zeros(nx,1)];
               
        %Jecon  = -Uk(1)*(2*Xk(2) - 0.5)* delta_time;
        %Jecon  = -Uk(1)* delta_time;
        %Jstate =(Qmax(1:nx,1).*(Xk - xdot_val_rf_ss))' * (Xk - xdot_val_rf_ss) * delta_time;
%         if check == 0
%             Jstate = 2*(Xk - xs(1:3))' * (Xk - xs(1:3)) * delta_time;
%         else
%             Jstate = 2*(Xk - xf)' * (Xk - xf) * delta_time;
%         end     
        if check == 0
            %Jstate = (Xk - xs(1:3))' * (Xk - xs(1:3)) * delta_time;
            Jstate = 0;
        else
            %Jstate = (Qc(1:3,1).*(Xk - xf))' * (Xk - xf) * delta_time;
            Jstate = (Xk - xf)' * (Xk - xf) * delta_time;
        end  

        Jecon  = -Xk(2)* delta_time;

        
%         % compute rotated cost function
%         fm  = f(Xk,Uk);
%         % load Lagrange multipliers from steady-state optimization
%         %Jmodel = le'*fm;
%         Jmodel = le'*(Xk - fm);
        
        %J = J + Jecon;
        J = J + Jecon + Jcontrol;
        %J = J + Jecon + Jcontrol + Jstate;

    end
end