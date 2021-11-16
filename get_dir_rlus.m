function [pi_sol] = get_dir_rlus(c,A_eq,r)
    options = optimoptions('linprog','Display','none');
	pi_sol = linprog(c,[],[],A_eq,ones(2*r,1),zeros(r*r,1),ones(r*r,1),options);
    pi_sol = reshape(pi_sol,[r,r])';
end
