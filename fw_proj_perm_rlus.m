function [pi_sol] = fw_proj_perm_rlus(P_hat,A_eq)
    options = optimoptions('linprog','Display','none');
    n      = size(P_hat,1);
    c      = reshape(P_hat,[n^2,1]);
	pi_sol = linprog(-c,[],[],A_eq,ones(2*n,1),zeros(n*n,1),ones(n*n,1),options);
    pi_sol = reshape(pi_sol,[n,n])';
end