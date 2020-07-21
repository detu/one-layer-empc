function [Tall, xmeasureAll, uAll, xuSSAll, primalNLP, params, runtime] = iNmpcCyclic(optProblem, system, mpciterations, N, T, tmeasure, xmeasure, u0, varargin)
%INMPCCYCLIC Summary of this function goes here
% 
% [OUTPUTARGS] = INMPCCYCLIC(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2018/01/17 21:39:15 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2018

import casadi.*

Tall    = [];
Xall    = zeros(mpciterations,size(xmeasure,1));
Uall    = zeros(mpciterations,size(u0,1));

xmeasureAll = [];
uAll        = [];
runtime     = [];
x_nlp_opt   = [];
u_nlp_opt   = [];
mpciter     = 1;
xuSSAll     = [];
xuSSt        = [xmeasure;u0(1)];
%load noise1pct.mat;

global nx delta xf uf check le;
xf = xmeasure;
uf = u0(1);

while(mpciter <= mpciterations)
    
    fprintf('-----------------------------\n');
    fprintf('MPC iteration: %d\n', mpciter);
    
    % change parameter value
    if mpciter == 1
        check = 0;
        delta = 0.55;
    end
    
    if mpciter == 25
        check = 0;
        delta = 1.1;
    end
    
    if mpciter == 80
        check = 0;
        delta = 2.2;
    end 
    
    % obtain new initial value
    [t0, x0] = measureInitialValue ( tmeasure, xmeasure );
    
    x0 = max(min(x0,1.0),0); % restrict to boundaries
    
    % first trial is without measurement noise
    x0_measure         =  x0;    % without noise
    
    % check constraints on boundary
    %x0_measure = max(min(x0_measure,1.0),0); % restrict to boundaries
     
    % ideal NMPC:
    [primalNLP, dualNLP, lb, ub, ~, params, elapsedtime] = solveOptimalControlProblem(optProblem, system, N, t0, x0, u0, T, mpciter, u_nlp_opt, x_nlp_opt, x0_measure, xuSSt);
    
    % re-arrange NLP results
    %[u_nlp_opt, x_nlp_opt] = plotStatesN(primalNLP, lb, ub, N);
    [u_nlp_opt, x_nlp_opt, xuSS] = arrangeOptResults(primalNLP, lb, ub, N);
    %xuSSt = xuSS';
    xuSSt = [x_nlp_opt(end-2);x_nlp_opt(end-1);x_nlp_opt(end);u_nlp_opt(end)];
%     if mpciter == 1
%         check = 1;
%         xf = [xuSS(1); xuSS(2); xuSS(3)];
%         uf = xuSS(4);
%         %save xuf.mat xf uf;
%     end
    
    if mpciter == 1 || mpciter == 25 || mpciter == 80
    %if mpciter == 1
        check = 1;
        xf = [xuSS(1); xuSS(2); xuSS(3)];
        uf = xuSS(4);
        %save xuf.mat xf uf;
        
%         % Compute Greshgoring bound here !
%         %% Compute Hessian and perform Greshgorin convexification
%         u1      = MX.sym('u1');     % u  - heat flux through cooling jacket
%         x1      = MX.sym('x1');     % x1 - concentration of R
%         x2      = MX.sym('x2');     % x2 - concentration of P1
%         x3      = MX.sym('x3');     % x3 - concentration of P2
%         %epsilon = MX.sym('epsilon');
%         x       = [x1;x2;x3;u1];
%         x1dot   = 1 - 1e4*x1^2*exp(-1/x3) - 400*x1*exp(-delta/x3) - x1;
%         x2dot   = 1e4*x1^2*exp(-1/x3) - x2;
%         x3dot   = u1 - x3;
%         ceq     = [x1dot;x2dot;x3dot];
%         % symbolic variable for dual variables
%         l1      = MX.sym('l1');
%         l2      = MX.sym('l2');
%         l3      = MX.sym('l3');
%         
%         xsol = primalNLP(1:4);
        lambda.eqnonlin = dualNLP(1:3);
        le   = lambda.eqnonlin;
%         % extract Lagrange multiplier
%         %l   = vertcat(l{:});
%         l  = [l1;l2;l3];
%         L  = -x2 + l'*ceq;
%         
%         Lagr = Function('Lagr', {x,l}, {L}, char('x','l'), char('Lagr'));
%         %Jac  = Function(Lagr.jacobian('x','Lagr'));
%         H    = Function(Lagr.hessian('x','Lagr'));
%         cons = Function('Const', {x}, {ceq}, char('x'), char('cons'));
%         Jcon = Function(cons.jacobian('x','cons'));
%         
%         eqVal = cons(xsol);
%         eqVal = full(eqVal);
%         Hx   = H(xsol,lambda.eqnonlin);
%         Hx   = full(Hx);
% %         Jac  = Jcon(xsol);
% %         Jac  = full(Jac);
% %         rH   = null(Jac)'*Hx*null(Jac)
% %         erH  = eig(rH)
%         
%         
%         [Hxxl,Qc]   = Greshgorin(Hx);
%         %save QmaxCstr.mat Qmax;
    end
    
    
    % save open loop solution for error computation
    %iNmpcData(mpciter).z = x_nlp_opt;

    z1 = x_nlp_opt(1:nx,5);  % 5 = (d+1) + 1 (d=number of collocation point)
    
    
    % record information
    Tall                = [Tall t0];
    Xall(mpciter+1,:)   = x0';
    Uall(mpciter+1,:)   = u0(:,1);
    
    
    % Apply control to process with optimized control from path-following
    % algorithm. 
    %x0 = xmeasure;  % from the online step 
    x0 = x0_measure; % NEW CHANGE 28.09.2017
    [tmeasure, xmeasure] = applyControl(system, T, t0, x0, u_nlp_opt);
    
%     %ObjVal(mpciter) = computeObjectiveFunctionValues(u_nlp_opt(:,1),xmeasure); % 
%     %ObjVal(mpciter) = mpciter; % DUMMY!
%     ObjVal(mpciter) = computeObjFuncCstr(u_nlp_opt(:,1),xmeasure); % CSTR only
    
    % store output variables
    xmeasureAll = [xmeasureAll;xmeasure];
    uAll        = [uAll;u_nlp_opt(:,1)];   
    runtime     = [runtime;elapsedtime];
    xuSSAll     = [xuSSAll;xuSS];
    
    % prepare restart
    u0 = shiftHorizon(u_nlp_opt);
    
    mpciter = mpciter+1;
end
xmeasureAll = reshape(xmeasureAll,nx,mpciterations);    

end

function [t0, x0] = measureInitialValue ( tmeasure, xmeasure )
    t0 = tmeasure;
    x0 = xmeasure;
end

function [tapplied, xapplied] = applyControl(system, T, t0, x0, u0)
    xapplied = dynamic(system, T, t0, x0, u0(:,1));
    %xapplied = dynamic(system, T, t0, x0, round(u0(:,1),1));
    tapplied = t0+T;
end

function u0 = shiftHorizon(u)
    u0 = [u(:,2:size(u,2)) u(:,size(u,2))];
end

function [u, lamda, lbw, ubw, objVal, params, elapsednlp] = solveOptimalControlProblem(optProblem, system, N, t0, x0, u0, T, mpciter, u_nlp, x_nlp, x0_measure, xuSS)

    import casadi.*
    % call ODE15s N-times initial guess in optimization
    x(1,:) = x0';
    for k=1:N
        x(k+1,:) = x0';    % ONE-TIME SIMULATION !
    end

    [J,g,w0,w,lbg,ubg,lbw,ubw,params] = optProblem(x, u0, N, x0_measure, xuSS);
    prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
    options = struct;
    options.ipopt.tol                = 1e-12;
    options.ipopt.constr_viol_tol    = 1e-12;  
    solver = nlpsol('solver', 'ipopt', prob, options);
    
    % Solve the NLP
    startnlp   = tic;
    sol        = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
    elapsednlp = toc(startnlp);
    fprintf('IPOPT solver runtime = %f\n',elapsednlp);
    success = strcmp(solver.stats.return_status,'Infeasible_Problem_Detected');
    if (success)
        keyboard;
    end

    u      = full(sol.x);
    lamda  = full(sol.lam_g);
    objVal = full(sol.f);
end

function [x, t_intermediate, x_intermediate] = dynamic(system, T, t0, x0, u0)
    x = system(t0, x0, u0, T);
    x_intermediate = [x0; x];
    t_intermediate = [t0, t0+T];
end

function [H,Q] = Greshgorin(H)
numH    = size(H,1);
Q       = zeros(numH,numH);
%Q       = eye(numH);
%delta   = 1e-1;
delta   = 1e-2;
%delta   = 2.5e-1;  % normal case
%delta   = 0.5;
%delta   = 1;
%delta   = 1.5;
%delta   = 2;
%delta   = 2.5;      % with measurement noise 1 percent
%delta   = 5;
%delta   = 100;
%delta   = 1e-6;
for i=1:numH  % iterate all row of Hessian
    sumRow = 0;
    for j=1:numH
        if j ~= i
            sumRow = sumRow + abs(H(i,j));
        end
    end
    %sumRow = sumRow - abs(H(i,i));
    
    if H(i,i) <= sumRow   % include equality 
    %if H(i,i) < sumRow 
        %Q(i,i) = sumRow - H(i,i) + delta;
        Q(i,i) = sumRow - H(i,i) + delta;
    end
end

% % loop over Qm to obtain maximum number
% for i=1:numH
%     if Q(i,i) > Qm(i,1)
%         Qm(i,1) = Q(i,i);
%     end
% end
Q = diag(Q);
end
