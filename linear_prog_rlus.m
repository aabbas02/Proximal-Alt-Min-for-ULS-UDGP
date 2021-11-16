function [pi_hat] = linear_prog_rlus(r,Y_hat,Y)
    options = optimoptions('linprog','Display','none');
    n      = size(Y_hat,1);
    pi_hat = zeros(n);
    A_eq = zeros(2*r,r*r);
    for i = 1 : r
        A_eq(i,(i-1)*r+1:i*r) = 1;
    end
    for i = 1 : r
        A_eq(i+r,i:r:i+(r-1)*r) = 1;
    end
    b_eq = ones(2*r,1);
    for i = 1 : n/r
        c      = -reshape(Y((i-1)*r+1:i*r,:)*Y_hat((i-1)*r+1:i*r,:)',[1,r^2]);        
        pi_sol = linprog(c,[],[],A_eq,b_eq,zeros(r*r,1),ones(r*r,1),options);
        pi_hat((i-1)*r+1:i*r,(i-1)*r+1:i*r) = reshape(pi_sol,[r,r])';
    end   
end
