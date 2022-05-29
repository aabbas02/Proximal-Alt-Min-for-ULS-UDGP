str = pwd;
k   = strfind(str,'\');
str = str(1:k(end));
cd (str)
addpath('datasets','code\prox_alt_min','code\misc','code\benchmarks\rlus')
clc
close all;
MC              = 1;
m               = 1;
n_              = [2016];
r_              = [96];
lbd             = 100;
d               = 10;
[ha, ~] = tight_subplot(1,3,[.01 .01],[.1 .01],[.01 .01])
for j = 1 : length(r_)
	r           = r_(j);
    for i = 1 : length(n_)
        n          = n_(i);
        for k = 1 : MC
               [Y] = get_pca_col_yale(d); %each column is data point
               [B,~,~]=svd(Y,'econ');
			   B = B(:,1:d);
               Y   = Y(:,16); % draw random data point, in this case 16th column 
               pi_ = make_r_local_permutation(n,r);
               axes(ha(1));
               plot_col_yale(Y)            
               Y_permuted      = pi_*Y;
               %title('Scrambled')
               axes(ha(2));
               plot_col_yale(Y_permuted)               
               [energy,pi_lp]    = lp_ls_alt_min_prox(B,Y_permuted,r,lbd);
               d_H              = sum(sum(pi_ ~= pi_lp))/2;
               %title('Reconstructed')
               axes(ha(3));
               plot_col_yale(pi_lp'*Y_permuted)
        end
    end
end
p_SNR = 10*log10(1/(norm(Y/max(Y) - pi_lp'*Y_permuted/max(Y_permuted),'fro')^2/n))
fig = gcf;
exportgraphics(fig,'yale.pdf','Resolution',300) 
function [data,B] = get_pca_col_yale(r)
    [data] = load_yale_compressed;
    data = data(:,1:64);
    [B,~,~]=svd(data,'econ');
    B = B(:,1:r);
end
function plot_col_yale(img)
   % imshow(reshape(round(img),[192,168]),[0,255]);
   imshow(reshape(round(img),[48,42]),[0,255]);
end
function [x] = load_yale_compressed
    load('yale_compressed.mat')
    x = zeros(2016,64);
    x_ = zeros(48,42);
    for i = 1 : 128 
       x_(1:48,1:42) = faces(i,1:2:96,1:2:84); %down sample by 2x2
       x(:,i) = reshape(x_,2016,1);
    end 
end


