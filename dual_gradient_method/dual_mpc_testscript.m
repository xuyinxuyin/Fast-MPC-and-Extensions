
%% Parameters
T  = 10;
n = 12; % state dimension
m = 3; % input dimension
l = 18;
k = 24;

% objective matrices
Q = eye(n);
R = eye(m);
Qf = Q;
S  = [];
q  = [];
r  = [];
qf = [];

% dynamics
k = 1;          % spring constant
lam = 0;        % damping constant
Aa = -2*k;
Ab = -2*lam;
Ac = k;
Ad = lam;

Acts = [zeros(n/2) eye(n/2);
        [Aa,Ac,0,0,0,0,Ab,Ad,0,0,0,0;
         Ac,Aa,Ac,0,0,0,Ad,Ab,Ad,0,0,0;
         0,Ac,Aa,Ac,0,0,0,Ad,Ab,Ad,0,0;
         0,0,Ac,Aa,Ac,0,0,0,Ad,Ab,Ad,0;
         0,0,0,Ac,Aa,Ac,0,0,0,Ad,Ab,Ad;
         0,0,0,0,Ac,Aa,0,0,0,0,Ad,Ab]];

Bcts = [zeros(n/2,m);
        [1, 0, 0;
        -1, 0, 0;
         0, 1, 0;
         0, 0, 1;
         0,-1, 0;
         0, 0,-1]];

% convert to discrete-time system
% sampling time
ts = 0.5;
A = expm(ts*Acts);
B = (Acts\(expm(ts*Acts)-eye(n)))*Bcts;

% physical limits
Xmax = 4;
Umax = 0.5;
xmin = -Xmax*ones(n,1);
xmax = Xmax*ones(n,1);
umin = -Umax*ones(m,1);
umax = Umax*ones(m,1);

w = 2*rand(n,1)-1;
w(1:n/2,:) = 0;
w = 0.5*w;

% initial state
x0 = zeros(n,1);
% final state
xf = [];
kappa = 0.01;
Kmax = 5;

test = Fast_MPC(Q,R,S,Qf,q,r,qf,xmin,xmax,umin,umax,T,x0, A,B,w,xf, []);   % Build class
[H, g] = form_objective_function(test);
[P, h] = form_inequality_const(test);
[C, b] = form_equality_const(test);

% Coefficients of dual problem
H_inv = inv(H);
dual_A = [P*H_inv*(P'), P*H_inv*(C');
      C*H_inv*(P'), C*H_inv*(C')]/4;
dual_B = [P*g/2 + h; C*g/2 + b];
dual_C = g'*H_inv*g/4;

%% Solving
fprintf('\n\nSingle MPC Step Computation Time Comparision\n');
fprintf('======================================================================\n');
% Solving with infeasible newton method, solve exact problem in each MPC step
tic;
fprintf('Solving with infeasible newton method -- full version\n');
[x_opt_full] = test.mpc_solve_full;
t_full = toc;
fprintf('Infeasible start newton exact problem finished in %3f sec\n',t_full);

% Fixed log barrier + fixed newton step
tic;
fprintf('Solving with infeasible newton method -- fixed kappa and max newton steps\n');
[x_opt_lgnw] = test.mpc_fixed_log_newton(Kmax,kappa);
t_both = toc;
fprintf('Infeasible start newton with both fixed kappa and Kmax finished in %3f sec\n',t_both);

% Fixed log barrier + fixed gradient step
tic;
fprintf('Solving with gradienet method -- fixed kappa and max gradient steps\n');
[x_grad, u_grad, ~] = dual_mpc(dual_A, dual_B, dual_C, kappa, T, n, m, H_inv, g, P, C);
t_grad = toc;
fprintf('dual gradient method with both fixed kappa and Kmax finished in %3f sec\n',t_both);

fprintf('======================================================================\n');

%% Plotting
x_full = zeros(T*n,1);
u_full = zeros(T*m,1);
for i=1:(m+n):length(x_opt_full)
    if i==1
        u_full(i:i+m-1) = x_opt_full(i:i+m-1);
        x_full(i:i+n-1) = x_opt_full(i+m:i+m+n-1);
    else
        u_full((i-1)/(m+n)*m+1:(i-1)/(m+n)*m+m) = x_opt_full(i:i+m-1);
        x_full((i-1)/(m+n)*n+1:(i-1)/(m+n)*n+n) = x_opt_full(i+m:i+m+n-1);
    end
end

x_lgnw = zeros(T*n,1);
u_lgnw = zeros(T*m,1);
for i=1:(m+n):length(x_opt_lgnw)
    if i==1
        u_lgnw(i:i+m-1) = x_opt_lgnw(i:i+m-1);
        x_lgnw(i:i+n-1) = x_opt_lgnw(i+m:i+m+n-1);
    else
        u_lgnw((i-1)/(m+n)*m+1:(i-1)/(m+n)*m+m) = x_opt_lgnw(i:i+m-1);
        x_lgnw((i-1)/(m+n)*n+1:(i-1)/(m+n)*n+n) = x_opt_lgnw(i+m:i+m+n-1);
    end
end

figure();
subplot(2,1,1);
stairs(x_full); hold on;
stairs(x_lgnw);
stairs(x_grad);
legend(...
       ['exact newton (' num2str(t_full * 1000, '%.2f') 'ms)'],...
       ['fixed kappa + Kmax (' num2str(t_both * 1000, '%.2f') 'ms)'], ...
       ['fixed kappa + Gmax (' num2str(t_grad * 1000, '%.2f') 'ms)']);
xlabel('t');
ylabel('$x_1(t)$', 'Interpreter', 'latex');
axis tight;
title('Predicted $x_1(t)$ in the first MPC iteration', 'Interpreter', 'latex');

subplot(2,1,2);
stairs(u_full); hold on;
stairs(u_lgnw);
stairs(u_grad);
legend(...
       ['exact newton (' num2str(t_full * 1000, '%.2f') 'ms)'],...
       ['fixed kappa + Kmax (' num2str(t_both * 1000, '%.2f') 'ms)'], ...
       ['fixed kappa + Gmax (' num2str(t_grad * 1000, '%.2f') 'ms)']);
xlabel('t');
ylabel('$u_1(t)$', 'Interpreter', 'latex');
axis tight;
title('Predicted $u_1(t)$ in the first MPC iteration', 'Interpreter', 'latex');


%% Supporting functions
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
    if (size(obj.x_min,1)~=size(obj.Q,1)) || (size(obj.x_max,1)~=size(obj.Q,1))
        error('Check the state inequality constraints dimensions');
    end
    if (size(obj.u_min,1)~=size(obj.R,1)) || (size(obj.u_max,1)~=size(obj.R,1))
        error('Check cotrol iequality constraint dimension');
    end

    %% Parameters
    x_min = obj.x_min;
    x_max = obj.x_max;
    u_min = obj.u_min;
    u_max = obj.u_max;
    T = obj.T;
    n = size(x_min,1);
    m = size(u_min,1);
    x = obj.x0;

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
    h(i:2*(m+n)+i-1) = [u_max;-u_min;x_max;-x_min];
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
    n = size(A,2);
    m = size(B,2);
    x = obj.x0;
    u = obj.R;
    T = obj.T;
    C = zeros(T*n,T*(n+m));
    b = zeros(T*n,1);
    xf = obj.x_final;

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
