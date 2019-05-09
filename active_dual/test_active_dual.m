
%% Parameters

 %n = 8;                      % Dimension of state
 %m = 3;                      % Dimension of control
n=3;
m=2;
Q = eye(n);                 % State stage cost
R = eye(m);                 % Control stage cost
S = [];                     % State control coupled cost
%Qf = 50*eye(n);             % Terminal state cost
q = [];                     % Linear state cost
r = [];                     % Linear control cost
qf = [];                    % Terminal state cost
Xmax =1;                   % State upper limit
Umax =1;                 % Control upper limit
xmin = -Xmax*ones(n,1);     % State lower bound
xmax = Xmax*ones(n,1);      % State upper bound
umin = -Umax*ones(m,1);     % Cotrol lower bound
umax = Umax*ones(m,1);      % Control upper bound



high_limit = 1;
low_limit = 0;
%A = (high_limit-low_limit).*rand(n,n) + ones(n,n)*low_limit;    % Random A (State transition) matrix
%B = (high_limit-low_limit).*rand(n,m) + ones(n,m)*low_limit;    % Random B (Control matrix) matrix
A=load('A.mat');
A=A.A;
B=load('B.mat');
B=B.B;
% A=0.5*ones(n,n);
% B=0.5*ones(n,m);
A = A./(max(abs(eig(A))));      % Spectral radius of A within 1

high_limit_w = 1;
low_limit_w = 0;
w = 0*(high_limit_w-low_limit_w).*rand(n,1) + ones(n,1)*low_limit_w;  % Random noise vector

T = 30;                         % Horizon length
%T=3;
% xf = 1*ones(n,1);                 % Initial state (random)
% x0 = rand(n,1);              % Terminal state
x0=1*ones(n,1);
%xf=rand(n,1);
xf=load('xf.mat');
xf=xf.xf;
[x_mat,u_mat]=active_dual(Q,R,xmin,xmax,umin,umax,T,x0,A,B,w,xf);


% x_mat = zeros(T*n,1);
% u_mat = zeros(T*m,1);
% for i=1:(m+n):length(x_opt_mat)
%     if i==1
%         u_mat(i:i+m-1) = x_opt_mat(i:i+m-1);
%         x_mat(i:i+n-1) = x_opt_mat(i+m:i+m+n-1);
%     else
%         u_mat((i-1)/(m+n)*m+1:(i-1)/(m+n)*m+m) = x_opt_mat(i:i+m-1);
%         x_mat((i-1)/(m+n)*n+1:(i-1)/(m+n)*n+n) = x_opt_mat(i+m:i+m+n-1);
%     end
% end
x_mat=x_mat(1:n:(T-1)*n+1);
u_mat=u_mat(1:m:(T-1)*m+1);
figure(1);
stairs(x_mat)
ylabel('x_t(1)')
title('dual problem for MPC solved with active set')
figure(2);
stairs(u_mat);
ylabel('u_t(1)')
title('dual problem for MPC solved with active set')

