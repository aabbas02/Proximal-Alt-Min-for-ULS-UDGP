%This script generates fig. 3(c) 
str = pwd;
k   = strfind(str,'\');
str = str(1:k(end));
cd (str)
addpath('code\misc','code\prox_alt_min') 
clc
close all;
clear all
MC              = 25;
SNR             = 100;
d               = 100;
r               = 100;
m_              = [10 15 20 25 30 35 40 45 50];
n_              = [1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500];
err_X           = zeros(length(n_),length(m_));
err_Pi          = zeros(length(n_),length(m_));
options         = optimoptions('linprog','Display','none');
lbd             = 100;
X_full          = randn(d,max(m_));
for i = 1 : length(n_)
    n   = n_(i);
    pi_ = make_r_local_permutation(n,r);
for j = 1 : length(m_)
    m           = m_(j);
    X           = X_full(:,1:m);
    for k = 1 : MC
        B                    = randn(n,d);
        Y                    = B*X;  
        noise_var            = norm(X,'fro')^2  / (SNR*m);
        W                    = sqrt(noise_var)*randn(n,m);
        Y_permuted           = pi_*Y;
        Y_permuted_noisy     = Y_permuted + W;
        [energy,pi_lp]       = lp_ls_alt_min_prox(B,Y_permuted_noisy,r,lbd);
        d_H                  = sum(sum(pi_ ~= pi_lp))/(2*n);
        err_Pi(i,j)          = d_H + err_Pi(i,j);
    end
end
end
err_Pi = err_Pi/MC;
imagesc(flip(err_Pi))
colormap(gray)
set(gca,'ytick',1:length(n_))
set(gca,'xtick',1:length(n_))
set(gca,'yticklabels',flip(n_))
set(gca,'xticklabels',m_)
xlabel('$m$','Interpreter','Latex','Fontsize',15)
ylabel('$n$','Interpreter','Latex','Fontsize',15)
colorbar
title('Fraction Hamming distortion $d_H/n$','Interpreter','Latex','Fontsize',15) 
