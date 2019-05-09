classdef Fast_MPC
% FAST_MPC Problem formulation
% min sum_{t0}^{t0+T} [x(t) u(t)'] [Q S; S' R] [x(t) u(t)] + q'x(t) + r'u(t)
% s.t.
%     x(t+1) = Ax(t) + Bu(t) + w(t), t = t0, t0+1, ..., t0+T
%     x(t) <= xmax, t = t0, t0+1, ..., t0+T
%     x(t) >= xmin, t = t0, t0+1, ..., t0+T
%     u(t) <= umax, t = t0, t0+1, ..., t0+T
%     u(t) >= umin, t = t0, t0+1, ..., t0+T
%     x(t0) = x0
%     x(t0+T) = xf
%
% Note x(t0+T) = xf will be encoded into the equality constrait and  x(t0) = x0
% will be implictly satisfied by choosing x0 as the starting point.
properties
    T;
    Q; R; S; q; r; Qf; qf;
    A; B; w;
    xmin; xmax; umin; umax;
    x0; xf;
    n; m;
end
methods
    function obj = ...
        Fast_MPC(T, Q, R, S, q, r, Qf, qf, A, B, w, xmin, xmax, umin, umax, x0, xf)
        obj.T = T;

        obj.Q = Q;
        obj.R = R;
        obj.S = S;
        obj.q = q;
        obj.r = r;
<<<<<<< HEAD
        obj.qf = qf; 
        obj.x_min = xmin;
        obj.x_max = xmax;
        obj.u_min = umin;
        obj.u_max = umax;
        obj.T = T;
        obj.x0 = x0;
=======
        obj.Qf = Qf;
        obj.qf = qf;

>>>>>>> fc9a6298c8896c2d8a2da90a38feb5c19ff3dd5a
        obj.A = A;
        obj.B = B;
        obj.w = w;

        obj.xmin = xmin;
        obj.xmax = xmax;
        obj.umin = umin;
        obj.umax = umax;

        obj.x0 = x0;
        obj.xf = xf;

        obj.n = length(xmin);
        obj.m = length(umin);
    end

    % matlab fmincon, solve exact problem in each MPC step
    function [x_opt] = matlab_solve(obj)
        kappa = 1;
        mu = 1/10;
        options = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');
        z0 = zeros(obj.T * (obj.n + obj.m),1);
        while kappa*length(z0) >= 10e-3
            % formulate MPC
            [P,h] = form_inequality_const(obj);
            [H,g] = form_objective_function(obj);
            J = @(z)(z'*H*z + g'*z + kappa*(-sum(log(h - P*z))));
            [A_eq,b_eq] = form_equality_const(obj);

            x_opt = fmincon(J,z0,[],[],A_eq,b_eq,[],[],[],options);
            z0 = x_opt;
            kappa = mu*kappa;
        end
    end

    % infeasible newton method, solve exact problem in each MPC step
    function [x_opt] = mpc_solve_full(obj)
        kappa = 1;
        mu = 1/10;
        [H,g] = form_objective_function(obj);
        [P,h] = form_inequality_const(obj);
        [C,b] = form_equality_const(obj);
        z = zeros(size(g));
        while kappa*length(z) >= 10e-3
            x_opt = infeasible_newton_solver(H,g,P,h,C,b,kappa,z,[]);
            kappa = mu*kappa;
            z = x_opt;
        end

    end

    % fixed barrier parameter
    function [x_opt] = mpc_fixed_log(obj,kappa)
        [H,g] = form_objective_function(obj);
        [P,h] = form_inequality_const(obj);
        [C,b] = form_equality_const(obj);
        z = zeros(size(g));
        x_opt = infeasible_newton_solver(H,g,P,h,C,b,kappa,z,[]);
    end

    % fixed iteration number
    function [x_opt] = mpc_fixed_newton(obj,max_nt_iter)
        kappa = 1;
        mu = 1/10;
        [H,g] = form_objective_function(obj);
        [P,h] = form_inequality_const(obj);
        [C,b] = form_equality_const(obj);
        z = zeros(size(g));
        while kappa*length(z) >= 10e-3
            x_opt = infeasible_newton_solver(H,g,P,h,C,b,kappa,z,max_nt_iter);
            kappa = mu*kappa;
            z = x_opt;
        end
    end

    function [x_opt] = mpc_fixed_log_newton(obj,max_nt_iter,kappa)
        [H,g] = form_objective_function(obj);
        [P,h] = form_inequality_const(obj);
        [C,b] = form_equality_const(obj);
        z = zeros(size(g));
        x_opt = infeasible_newton_solver(H,g,P,h,C,b,kappa,z,max_nt_iter);
    end

end
end

%% Supporting fuctions

function [H,g] = form_objective_function(obj)
%FAST_MPC_OBJECTIVE Return the objective function in matrix form
%    [H,g] = fast_mpc_objective(obj)
%    where the objective function is f = z'Hz + g'z

    %% Parameters
    T = obj.T;                 % Horizon time
    n = size(obj.Q,1);         % State vector dimension
    m = size(obj.R,1);         % Control vector dimension
    z = zeros(T*(n+m),1);

    R  = obj.R;
    Q  = obj.Q;
    Qf = obj.Qf;
    q  = obj.q;
    r  = obj.r;
    qf = obj.qf;
    x  = obj.x0;

    %% Checking initial conditions
    if (size(obj.Q,1) ~= size(obj.Q,2)) ||...
            (size(obj.Qf,1) ~= size(obj.Qf,2))
        error('State stage cost must a square matrix');
    elseif size(obj.R,1) ~= size(obj.R,2)
        error('Control stage cost must a square matrix');
    end

    % Checking state and control linear stage cost size
    if isempty(obj.q) == 0
        if (size(obj.q,1) ~= size(obj.Q,1))
            error('Linear state cost needs to be a vector of size n');
        end
    else
        q = zeros(size(Q,1),1);
    end

    if isempty(obj.r) == 0
        if (size(obj.r,1) ~= size(obj.R,1))
            error('Linear control cost needs to be a vector of size n');
        end
    else
        r = zeros(size(R,1),1);
    end
    if isempty(obj.qf) == 0
        if (size(obj.qf,1) ~= size(obj.Q,1))
            error('State terminal linear cost needs to be a vector of size n');
        end
    else
        qf = zeros(size(Q,1),1);
    end

    %% Objective construction
    H = zeros(T*(n+m),T*(n+m));
    for i=m+1:(n+m):size(H,1)-n
        H(i:i+(n+m)-1,i:i+(n+m)-1) = [Q zeros(n,m); zeros(m,n) R];
    end
    H(1:m,1:m) = R;
    H(end-n+1:end,end-n+1:end) = Qf;

    g = zeros(T*(n+m),1);
    for i=m+1:(n+m):size(g,1)
        if i == T*(n+m) -n+1
            g(end-n+1:end) = qf;
        else
            g(i:i+(n+m)-1) = [q;r];
        end
    end
    g(1:m) = r;
end

function [P,h] = form_inequality_const(obj)
%FAST_MPC_INEQ_CONST Return the inequality constraints in matrix form
%    [P,h] = fast_mpc_ineq_const(obj)
%    where the inequality constrait is Pz <= h

    %% Initial condition check
    if (size(obj.xmin,1)~=size(obj.Q,1)) || (size(obj.xmax,1)~=size(obj.Q,1))
        error('Check the state inequality constraints dimensions');
    end
    if (size(obj.umin,1)~=size(obj.R,1)) || (size(obj.umax,1)~=size(obj.R,1))
        error('Check cotrol iequality constraint dimension');
    end

    %% Parameters
    T = obj.T;
    n = size(obj.xmin,1);
    m = size(obj.umin,1);

    %% Inequality constraint construction
    P = zeros(2*T*(n+m),T*(n+m));
    h = zeros(2*T*(n+m),1);

    for i=1:2*(m+n):size(P,1)-2*(m+n)+1
        if i==1
            P(i:i+2*(m+n)-1,i:i+(m+n)-1) = [eye(m) zeros(m,n);(-eye(m)) zeros(m,n);...
                zeros(n,m) eye(n); zeros(n,m) -(eye(n))];
        else
            P(i:i+2*(m+n)-1,(i+1)/2:(i+1)/2+(m+n)-1) = [eye(m) zeros(m,n);(-eye(m)) zeros(m,n);...
                zeros(n,m) eye(n); zeros(n,m) -(eye(n))];
        end
    end

    for i=1:2*(m+n):2*T*(m+n)-2*(m+n)+1
    h(i:2*(m+n)+i-1) = [obj.umax;-obj.umin;obj.xmax;-obj.xmin];
    end
end

function [C,b] = form_equality_const(obj)
%FAST_MPC_EQ_CONST Return the equality constraint in matrix form
%    [C,b] = fast_mpc_eq_const(obj)
%    where the equality constraint is Cz = b

    %% Parameters
    A = obj.A;
    B = obj.B;
    w = obj.w;
    n = obj.n;
    m = obj.m;
    x = obj.x0;
    u = obj.R;
    T = obj.T;
    C = zeros(T*n,T*(n+m));
    b = zeros(T*n,1);
    xf = obj.xf;

    %% Initial condition check
    if isempty(A)
        error('Define the state dynamics/equality constrained matrix');
    elseif isempty(B)
        error('Define the control dynamics/equality constrained matrix');
    end

    if size(A,2) ~= size(x,1)
        error('The equality state dynamics matrix size does not match');
    elseif size(B,2) ~= size(u,2)
        error('The equality control dynamics matrix size does not match');
    elseif isempty(w)
        w = zeros(n,1);
    end

    %% Equality constraint construction

    for i=n:n:T*n-n+1
        if i==n
            C(i+1:i+n,i+((i/n)-1)*(n+m+n):i+n+((i/n)-1)*(n+m+n)+m+n-1) = [-A -B eye(n)];
            b(i+1:i+n) = w;
        else
            C(i+1:i+n,((i/n)-1)*(n+m)+m+1:((i/n)-1)*(n+m)+m+m+n+n) = [-A -B eye(n)];
            b(i+1:i+n) = w;
        end
    end
    C(1:n,1:m+n) = [-B eye(n)];
    % C(end-n+1:end,end-(n+m+n)+1:end) = [zeros(n,n) zeros(n,m) eye(n)];
    b(1:n) = A*x + w;
    % b(end-n+1:end) = xf;

    % attach the final state requirement to the end, so that x = xf will also be
    % honored as an equality constraint.
    % if empty then we just ignore it.
    if isempty(xf)~=1
    b = [b;xf];
    C = [C;zeros(n,size(C,2))];
    end
    C(end-n+1:end,end-n+1:end) = eye(n);
end
