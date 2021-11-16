function [Pi_1] =  admm(B,Y,r,rho)
    n        = size(B,1);
    d        = size(B,2);
    m        = size(Y,2);
    A_eq = zeros(n,n);
    for i = 1 : n/r
        A_eq((i-1)*r+1:i*r,(i-1)*r+1:i*r) = ones(r,r);
    end
    A_eq = reshape(A_eq,[1,n^2]);
    B_bar = (1/d)*sum(B,2);
    Y_bar = (1/m)*sum(Y,2);
    P_B = B*inv(B'*B)*B';
    mu  = 0*ones(n,n);
    [fval1,Pi_1] = lp_r_admm(-Y_bar*B_bar',A_eq); % - sign solves maximization problem
    [fval2,Pi_2] = lp_r_admm(+Y_bar*B_bar',A_eq); % + sign solves minimization problems
    if(abs(fval2^2) > abs(fval1^2))
        Pi_1 = Pi_2;
    end
    Pi_2 = Pi_1;
    start = 1;
    while(norm(Pi_1 - Pi_2) ~= 0 || start==1)
        start = 0;
        C = -(Y*Y')*Pi_2*P_B' + mu - rho*Pi_2;
        [~,Pi_1] = lp_r_admm(C,A_eq);
        C = -(Y*Y')*Pi_1*P_B - mu - rho*Pi_1;
        [~,Pi_2] = lp_r_admm(C,A_eq);
        mu  = mu + rho*(Pi_1 - Pi_2);
        norm(Pi_2-Pi_1,'fro');
        %-trace(Pi_1*P_B*Pi_2'*(Y*Y'))
    end
end