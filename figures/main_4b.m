% download  MNIST "train set" from https://pjreddie.com/projects/mnist-in-csv/
addpath('datasets','code\prox_alt_min','code\misc')
clc
close all;
MC              = 1;
m               = 1;
n_              = [784];
r_              = [28];
lbd             = 100;
d               = 5;
[ha, ~] = tight_subplot(1,3,[.01 .01],[.1 .01],[.01 .01])
for j = 1 : length(r_)
	r           = r_(j);
    for i = 1 : length(n_) 
        n          = n_(i);
        for k = 1 : MC
               [Y] = get_mnist();
               [B,~,~]=svd(Y,'econ');
			   B = B(:,1:d);
               Y   = Y(:,125);
               pi_ = make_r_local_permutation(n,r);
               axes(ha(1));
               plot_img(Y)               
               Y_permuted      = pi_*Y;
               axes(ha(2));
               plot_img(Y_permuted)               
               [energy,pi_lp]    = lp_ls_alt_min_prox(B,Y_permuted,r,lbd);
               d_H              = sum(sum(pi_ ~= pi_lp))/2;
               axes(ha(3));
               plot_img(pi_lp'*Y_permuted)
        end
    end
end
p_SNR = 10*log10(1/ ( norm(Y/max(Y) - pi_lp'*Y_permuted/max(Y_permuted),'fro')^2/n))
fig = gcf;
exportgraphics(fig,'mnist.pdf','Resolution',300) 
function plot_img(img)
   % imshow(reshape(round(img),[192,168]),[0,255]);
   imshow(reshape(round(img),[28,28]),[0,255]);
end
function [Y] = get_mnist()
    load('mnist_train.csv')
    Y = mnist_train(:,2:785); % exclude labels
    Y = Y';
end


