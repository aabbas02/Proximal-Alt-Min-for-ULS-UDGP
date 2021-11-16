clc
close all;
MC              = 1;
m               = 1;
n               = (96*84/4);
r               = 96;
num_faces       = 10;
table           = zeros(2,64*num_faces);
options         = optimoptions('linprog','Display','none');
lbd             = 100;
d               = 10;
num_eg          = 64;
[Y_all] = get_pca_col_yale(num_faces); %each column is data point
pi_     = make_r_local_permutation(n,r);
t = 1;
	for i = 1 : num_faces
		Y = Y_all(:,(i-1)*num_eg+1:i*num_eg);
		[B,~,~]=svd(Y,'econ');
		B = B(:,1:d);
		for j = 1 : 64
			y   = Y(:,j);
			y_permuted        = pi_*y;
			% proposed
			[~,pi_lp]    	  = lp_ls_alt_min_prox(B,y_permuted,r,lbd);
            energy            = 10*log10(1/(norm(y/max(y) - pi_lp'*y_permuted/max(y_permuted),'fro')^2/n));
			table(1,t) 		  = energy;
            sum(table(1,:))
			% rlus
			[pi_fw]        	  = rlus(B,y_permuted,r);
			energy            = 10*log10(1/ ( norm(y/max(y) - (pi_fw*y_permuted)/max(y_permuted),'fro')^2/n));
			table(2,t) 		  = energy;
            sum(table(2,:))
            t = t+1;
		end
		i
	end
%save('yale_psnr_std_d_5.mat');
std_1 = std(table(1,:))
std_2 = std(table(2,:))
table = sum(table,2);
table = table/(num_eg*num_faces)
function [data] = get_pca_col_yale(num_faces)
    [data] = load_yale_compressed(num_faces);
    data = data(:,1:num_faces*64);
end
function [x] = load_yale_compressed(num_faces)
    load('yale_compressed.mat')
    x = zeros(2016,num_faces*64);
    x_ = zeros(48,42);
    for i = 1 : num_faces*64
       x_(1:48,1:42) = faces(i,1:2:96,1:2:84);
       x(:,i) = reshape(x_,48*42,1);
    end
end



