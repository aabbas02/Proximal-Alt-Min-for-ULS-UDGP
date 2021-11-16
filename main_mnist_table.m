% download  MNIST "train set" from https://pjreddie.com/projects/mnist-in-csv/
clc
close all;
clear all;
m               = 1;
n               = [784];
r               = [28];
num_digits      = 10;
options         = optimoptions('linprog','Display','none');
lbd             = 100;
d               = 5;
num_eg          = 200;
table           = zeros(2,num_digits*num_eg); % proposed,ojsp,admm
[Y_all,labels]  = get_pca_col(); 			  % each column is data point
pi_ = make_r_local_permutation(n,r);
t = 1;
for i = 0 : num_digits-1
		idx = find(labels == i);
		Y   = Y_all(:,idx);
		[B,~,~]=svd(Y,'econ');
		B = B(:,1:d);
        for j = 1 : num_eg
			y   = Y(:,j);
			y_permuted        = pi_*y;
			% proposed
			[~,pi_lp]    	  = lp_ls_alt_min_prox(B,y_permuted,r,lbd);
			energy            = 10*log10(1/(norm(y/max(y) - (pi_lp'*y_permuted)/max(y_permuted),'fro')^2/n));
			table(1,t) 		  =  energy;
            sum(table(1,:))
			% rlus
			[pi_fw]        	  = rlus(B,y_permuted,r);
			energy            = 10*log10(1/ ( norm(y/max(y) - (pi_fw*y_permuted)/max(y_permuted),'fro')^2/n));
			table(2,t) 		  = energy;
            sum(table(2,:))
			%-------------------------------------------------------------------------
            t = t+1;
        end
        i
end
%save('mnist_psnr_std.mat');
std_1 = std(table(1,:))
std_2 = std(table(2,:))
table = sum(table,2);
table = table/(num_eg*num_digits)
function [Y,labels] = get_pca_col()
    load('mnist_train.csv')
    Y = mnist_train(:,2:785); % exclude labels
    Y = Y';
    labels = mnist_train(:,1);
end



