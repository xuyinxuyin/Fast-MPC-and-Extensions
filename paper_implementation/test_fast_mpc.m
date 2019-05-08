% Testing FAST MPC class
clear all;
% close all;
%% Parameters

T = 10;                     % Horizon length
n = 8;                      % Dimension of state
m = 5;                      % Dimension of control
Q = eye(n);                 % State stage cost
R = eye(m);                 % Control stage cost
S = [];                     % State control coupled cost
Qf = 50*eye(n);             % Terminal state cost
q = [];                     % Linear state cost
r = [];                     % Linear control cost
qf = [];                    % Terminal state cost
Xmax = 10;                  % State upper limit
Umax = 2;                   % Control upper limit
xmin = -Xmax*ones(n,1);     % State lower bound
xmax = Xmax*ones(n,1);      % State upper bound
umin = -Umax*ones(m,1);     % Cotrol lower bound
umax = Umax*ones(m,1);      % Control upper bound

high_limit = 1;
low_limit = 0;
A = (high_limit-low_limit).*rand(n,n) + ones(n,n)*low_limit;    % Random A (State transition) matrix
B = (high_limit-low_limit).*rand(n,m) + ones(n,m)*low_limit;    % Random B (Control matrix) matrix

A = A./(max(abs(eig(A))));      % Spectral radius of A within 1

high_limit_w = 1;
low_limit_w = 0;
w = (high_limit_w-low_limit_w).*rand(n,1) + ones(n,1)*low_limit_w;  % Random noise vector

x0 = rand(n,1);                 % Initial state (random)
xf = 1*ones(n,1);               % Terminal state
test = Fast_MPC(Q,R,S,Qf,q,r,qf,xmin,xmax,umin,umax,T,x0, A,B,w,xf, []);   % Build class

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
stairs(x_mat);hold on;
stairs(x_full)
stairs(x_log);
stairs(x_nw);
stairs(x_lgnw);
legend(['fmincon (' num2str(t_mat, '%.2f') 's)'],...
       ['exact newton (' num2str(t_full * 1000, '%.2f') 'ms)'],...
       ['fixed kappa (' num2str(t_kappa * 1000, '%.2f') 'ms)'],...
       ['fixed Kmax (' num2str(t_kmax * 1000, '%.2f') 'ms)'],...
       ['fixed kappa + Kmax (' num2str(t_both * 1000, '%.2f') 'ms)']);
xlabel('t');
ylabel('$x_1(t)$', 'Interpreter', 'latex');
axis tight;
title('Predicted $x_1(t)$ in the first MPC iteration', 'Interpreter', 'latex');

subplot(2,1,2);
stairs(u_mat);hold on;
stairs(u_full)
stairs(u_log);
stairs(u_nw);
stairs(u_lgnw);
legend(['fmincon (' num2str(t_mat, '%.2f') 's)'],...
       ['exact newton (' num2str(t_full * 1000, '%.2f') 'ms)'],...
       ['fixed kappa (' num2str(t_kappa * 1000, '%.2f') 'ms)'],...
       ['fixed Kmax (' num2str(t_kmax * 1000, '%.2f') 'ms)'],...
       ['fixed kappa + Kmax (' num2str(t_both * 1000, '%.2f') 'ms)']);
xlabel('t');
ylabel('$u_1(t)$', 'Interpreter', 'latex');
axis tight;
title('Predicted $u_1(t)$ in the first MPC iteration', 'Interpreter', 'latex');
