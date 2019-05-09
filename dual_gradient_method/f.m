function val = f(x, dual_A, dual_B, dual_C, kappa, T, n, m)
if sum (x(1:2*T*(n+m))>0) == 2*T*(n+m)
    val = x'*dual_A*x + dual_B'*x + dual_C - kappa*sum(log(x(1:2*T*(n+m), 1)));
else
    val  = Inf;
end
end