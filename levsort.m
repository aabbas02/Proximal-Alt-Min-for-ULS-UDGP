function [pi_hat] = levsort(B,Y,r)    
    A_eq    = zeros(2*r,r*r);
    for i = 1 : r
        A_eq(i,(i-1)*r+1:i*r) = 1;
    end
    for i = 1 : r
        A_eq(i+r,i:r:i+(r-1)*r) = 1;
    end
	n           = size(Y,1);
    pi_hat  = zeros(n);
    options = optimoptions('linprog','Display','none');
    [Y,~,~] = svd(Y,'econ');
    [B,~,~] = svd(B,'econ');
    for i = 1:n/r
        %c    = reshape(diag(Y((i-1)*r+1:i*r,:)*Y((i-1)*r+1:i*r,:)')*...
        %               diag(B((i-1)*r+1:i*r,:)*B((i-1)*r+1:i*r,:)')'...
        %                     ,[1,r^2]);
        c     = reshape( diag(Y((i-1)*r+1:i*r,:)*Y((i-1)*r+1:i*r,:)')*...
                         diag(B((i-1)*r+1:i*r,:)*B((i-1)*r+1:i*r,:)')',...
                         [1,r^2]);
        %c    = reshape((Y((i-1)*r+1:i*r,:)*Y((i-1)*r+1:i*r,:)')*...
        %               (B((i-1)*r+1:i*r,:)*B((i-1)*r+1:i*r,:)')'...
        %                      ,[1,r^2]);
       temp = linprog(-c,[],[],A_eq,ones(2*r,1),zeros(r*r,1),[],options);
       pi_hat((i-1)*r+1:i*r,(i-1)*r+1:i*r) = reshape (temp,[r,r]);
    end
%     A_eq    = zeros(2*n,n*n);
%     for i = 1 : n
%         A_eq(i,(i-1)*n+1:i*n) = 1;
%     end
%     for i = 1 : n
%         A_eq(i+n,i:n:i+(n-1)*n) = 1;
%     end
%          c    = reshape(diag(Y*Y')*...
%                         diag(B*B')'...
%                         ,[1,n^2]);
%         temp = linprog(-c,[],[],A_eq,ones(2*n,1),zeros(n*n,1),[],options);
%         pi_hat = reshape (temp,[n,n]);
         
end