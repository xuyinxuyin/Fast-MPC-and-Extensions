function [X, U, T] = dual_mpc(dual_A, dual_B, dual_C, kappa, T, n, m, H_inv, g, P, C)

% parameters
epsilon = 1e-8;
alpha = 0.4;
beta = 0.8;

% dual variables
lambda = 1e-3*ones(2*T*(n+m), 1);
nu = 1e-3*ones(T*n, 1);
x = [lambda; nu];

% gradient method
x_grad = compute_gradient(x, dual_A, dual_B, T, n, m, kappa);
cnt = 0;
while(norm(x_grad)^2/2 > epsilon)
    if cnt > 5e2
        break
    end
    cnt = cnt + 1;
    
    % line search
    t = 1;
    x_next = x - t*x_grad;
    while(f(x_next, dual_A, dual_B, dual_C, kappa, T, n, m) > f(x, dual_A, dual_B, dual_C, kappa, T, n, m) - ...
            alpha*t*(x_grad'*x_grad))
        t = beta*t;
        x_next = x - t*x_grad;
    end
    
    x = x - t*x_grad;
    x_grad = compute_gradient(x, dual_A, dual_B, T, n, m, kappa);
end
z = -(2*H_inv)*(g + P'*x(1:2*T*(n+m), 1) + C'*x(end - T*n+1: end, 1));
z = reshape(z, [n+m, T]);

X = reshape(z(m+1: end, :), [T*n, 1]);
U = reshape(z(1:m, :), [T*m, 1]);
T = toc;
end