% Testing FAST MPC class

% Problem formulation
% min sum_{t0}^{t0+T} [x(t) u(t)'] [Q S; S' R] [x(t) u(t)] + q'x(t) + r'u(t)
% s.t.
%     x(t+1) = Ax(t) + Bu(t) + w(t), t = t0, t0+1, ..., t0+T
%     x(t) <= xmax, t = t0, t0+1, ..., t0+T
%     x(t) >= xmin, t = t0, t0+1, ..., t0+T
%     u(t) <= umax, t = t0, t0+1, ..., t0+T
%     u(t) >= umin, t = t0, t0+1, ..., t0+T
%     x(t0) = x0
%     x(t0+T) = xf

clear all;
% close all;
%% Parameters

T  = 10;
n = 12; % state dimension
m = 3; % input dimension

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

test = Fast_MPC(T, Q, R, S, q, r, Qf, qf, A, B, w, xmin, xmax, umin, umax, x0, xf);
%% Solving

fprintf('\n\nSingle MPC Step Computation Time Comparision\n');
fprintf('======================================================================\n');
% Native matlab solver
tic;
fprintf('Solving with matlab fmincon solver\n');
[x_opt_mat] = test.matlab_solve;
t_mat = toc;
fprintf('Matlab fmincon finished in %3f sec\n',t_mat);

% Solving with infeasible newton method, solve exact problem in each MPC step
tic;
fprintf('Solving with infeasible newton method -- full version\n');
[x_opt_full] = test.mpc_solve_full;
t_full = toc;
fprintf('Infeasible start newton exact problem finished in %3f sec\n',t_full);

% Fixed log barrier method k=0.01
tic;
kappa = 0.01;
fprintf('Solving with infeasible newton method -- fixed kappa\n');
[x_opt_log] = test.mpc_fixed_log(kappa);
t_kappa = toc;
fprintf('Infeasible start newton with fixed kappa = %f finished in %3f sec\n',kappa,t_kappa);

% Fixed newton step = 5
Kmax = 5;
tic;
fprintf('Solving with infeasible newton method -- fixed max newton steps\n');
[x_opt_nw] = test.mpc_fixed_newton(Kmax);
t_kmax = toc;
fprintf('Infeasible start newton with fixed K_max = %d finished in %3f sec\n',Kmax,t_kmax);

% Fixed log barrier + fixed newton step
tic;
fprintf('Solving with infeasible newton method -- fixed kappa and max newton steps\n');
[x_opt_lgnw] = test.mpc_fixed_log_newton(Kmax,kappa);
t_both = toc;
fprintf('Infeasible start newton with both fixed kappa and Kmax finished in %3f sec\n',t_both);

fprintf('======================================================================\n');

%% Plotting

x_mat = zeros(T*n,1);
u_mat = zeros(T*m,1);
for i=1:(m+n):length(x_opt_mat)
    if i==1
        u_mat(i:i+m-1) = x_opt_mat(i:i+m-1);
        x_mat(i:i+n-1) = x_opt_mat(i+m:i+m+n-1);
    else
        u_mat((i-1)/(m+n)*m+1:(i-1)/(m+n)*m+m) = x_opt_mat(i:i+m-1);
        x_mat((i-1)/(m+n)*n+1:(i-1)/(m+n)*n+n) = x_opt_mat(i+m:i+m+n-1);
    end
end

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

x_log = zeros(T*n,1);
u_log = zeros(T*m,1);
for i=1:(m+n):length(x_opt_log)
    if i==1
        u_log(i:i+m-1) = x_opt_log(i:i+m-1);
        x_log(i:i+n-1) = x_opt_log(i+m:i+m+n-1);
    else
        u_log((i-1)/(m+n)*m+1:(i-1)/(m+n)*m+m) = x_opt_log(i:i+m-1);
        x_log((i-1)/(m+n)*n+1:(i-1)/(m+n)*n+n) = x_opt_log(i+m:i+m+n-1);
    end
end

x_nw = zeros(T*n,1);
u_nw = zeros(T*m,1);
for i=1:(m+n):length(x_opt_nw)
    if i==1
        u_nw(i:i+m-1) = x_opt_nw(i:i+m-1);
        x_nw(i:i+n-1) = x_opt_nw(i+m:i+m+n-1);
    else
        u_nw((i-1)/(m+n)*m+1:(i-1)/(m+n)*m+m) = x_opt_nw(i:i+m-1);
        x_nw((i-1)/(m+n)*n+1:(i-1)/(m+n)*n+n) = x_opt_nw(i+m:i+m+n-1);
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
% stairs(x_mat);hold on;
stairs(x_full); hold on;
stairs(x_log);
stairs(x_nw);
stairs(x_lgnw);
% legend(['fmincon (' num2str(t_mat, '%.2f') 's)'],...
legend(...
       ['exact newton (' num2str(t_full * 1000, '%.2f') 'ms)'],...
       ['fixed kappa (' num2str(t_kappa * 1000, '%.2f') 'ms)'],...
       ['fixed Kmax (' num2str(t_kmax * 1000, '%.2f') 'ms)'],...
       ['fixed kappa + Kmax (' num2str(t_both * 1000, '%.2f') 'ms)']);
xlabel('t');
ylabel('$x_1(t)$', 'Interpreter', 'latex');
axis tight;
title('Predicted $x_1(t)$ in the first MPC iteration', 'Interpreter', 'latex');

subplot(2,1,2);
% stairs(u_mat);hold on;
stairs(u_full); hold on;
stairs(u_log);
stairs(u_nw);
stairs(u_lgnw);
% legend(['fmincon (' num2str(t_mat, '%.2f') 's)'],...
legend(...
       ['exact newton (' num2str(t_full * 1000, '%.2f') 'ms)'],...
       ['fixed kappa (' num2str(t_kappa * 1000, '%.2f') 'ms)'],...
       ['fixed Kmax (' num2str(t_kmax * 1000, '%.2f') 'ms)'],...
       ['fixed kappa + Kmax (' num2str(t_both * 1000, '%.2f') 'ms)']);
xlabel('t');
ylabel('$u_1(t)$', 'Interpreter', 'latex');
axis tight;
title('Predicted $u_1(t)$ in the first MPC iteration', 'Interpreter', 'latex');
