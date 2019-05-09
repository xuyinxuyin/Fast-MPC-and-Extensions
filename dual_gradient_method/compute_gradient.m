function grad = compute_gradient(x, dual_A, dual_B, T, n, m, kappa)
grad = 2*dual_A*x + dual_B - kappa*[1./x(1: 2*T*(n+m)); zeros(T*n, 1)];
end