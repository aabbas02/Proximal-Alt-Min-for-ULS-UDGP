function [pi_hat] = icml_20(B,Y,r)    
    A_eq    = zeros(2*r,r*r);
    for i = 1 : r
        A_eq(i,(i-1)*r+1:i*r) = 1;
    end
    for i = 1 : r
        A_eq(i+r,i:r:i+(r-1)*r) = 1;
    end
	n           = size(Y,1);
    pi_hat = zeros(n);
    options     = optimoptions('linprog','Display','none');
    for i = 1:n/r
        c           = reshape((Y((i-1)*r+1:i*r,:)*Y((i-1)*r+1:i*r,:)')*...
                              (B((i-1)*r+1:i*r,:)*B((i-1)*r+1:i*r,:)')...
                              ,[1,r^2]);
        temp = ...
        linprog(-c,[],[],A_eq,ones(2*r,1),zeros(r*r,1),[],options);
        pi_hat((i-1)*r+1:i*r,(i-1)*r+1:i*r) =...
        reshape(temp,[r,r]);
    end
end