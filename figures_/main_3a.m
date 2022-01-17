clc
close all;
clear all
str = pwd;
k   = strfind(str,'\');
str = str(1:k(end));
cd (str)
addpath('code\misc',...
        'code\prox_alt_min',...
        'code\benchmarks\rlus',...
        'code\benchmarks\biconvex',...
        'code\benchmarks\one_step',...
        'code\benchmarks\levsort') 

MC              = 5;
SNR             = 1000;
d               = 20;
m               = 20;
r_              = [2 5 10 20 25 40];
n               = 200;
d_H_levsort     = zeros(1,length(r_));
d_H_one_step    = zeros(1,length(r_));
d_H_rlus        = zeros(1,length(r_));
d_H_biconvex    = zeros(1,length(r_));
d_H_alt_min     = zeros(1,length(r_));
options         = optimoptions('linprog','Display','none');
rho_            = -3:1:5;
rho_            = 10.^rho_;
lbd             = 100;
for j = 1 : length(r_)
	r           = r_(j);
    for k = 1 : MC
                B                = randn(n,d);
                X                = randn(d,m);
                Y                = B*X;  
                noise_var   	 = norm(X,'fro')^2  / (SNR*m);
                W                = sqrt(noise_var)*randn(n,m);
                pi_              = make_r_local_permutation(n,r);
                Y_permuted       = pi_*Y;
                Y_permuted_noisy = Y_permuted + W;
				%---rlus  https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9440727 
                %t1_start = tic;
                [pi_fw]        = rlus(B,Y_permuted_noisy,r);
                %toc(t1_start)
                d_H            = sum(sum(pi_ ~= pi_fw))/(2*n);
                d_H_rlus(j)    = d_H + d_H_rlus(1,j);
				%---biconvex https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8849447
				d_H_min = 1;
                for i = 1 : length(rho_) % cross validate across rho paramter
                   rho              = rho_(i);
                   [pi_lp]          = admm(B,Y_permuted_noisy,r,rho);
                   d_H_             = sum(sum(pi_ ~= pi_lp))/(2*n);
                   if(d_H_ < d_H_min)
                      d_H_min = d_H_;
                   end
                end
                d_H_biconvex(j) = d_H_biconvex(j) + d_H_min;
                %---icml https://proceedings.mlr.press/v119/zhang20n.html
                pi_icml            = icml_20(B,Y_permuted_noisy,r);
                d_H_one_step(j)    = d_H_one_step(j) + sum(sum(pi_ ~= pi_icml))/(2*n);   
                %---levsort  https://people.eecs.berkeley.edu/~courtade/pdfs/DenoisingLinearModels_ISIT2017.pdf
                pi_lev             = levsort(B,Y_permuted_noisy,r);
                d_H_levsort(j)     = d_H_levsort(j) + sum(sum(pi_ ~= pi_lev))/(2*n);                   
				%---alt-min/proposed
                %t_alt_min = tic;
                [~,pi_lp]          = lp_ls_alt_min_prox(B,Y_permuted_noisy,r,lbd);
                %toc(t_alt_min)
                d_H                = sum(sum(pi_ ~= pi_lp))/(2*n);
                d_H_alt_min(j)     = d_H + d_H_alt_min(j); 
    end
    j
end
d_H_alt_min      = d_H_alt_min/MC;
d_H_one_step     = d_H_one_step/MC;
d_H_rlus         = d_H_rlus/MC;
d_H_levsort      = d_H_levsort/MC;
d_H_biconvex     = d_H_biconvex/MC;
styles =["r-x","g-x","b-x","k-x","m-x","c-x"];
hold on;
plot(r_,d_H_rlus,"g-x",...
    'DisplayName','RLUS',...
    'MarkerSize',11,'Linewidth',1.25);
plot(r_,d_H_biconvex,styles(3),'DisplayName','Biconvex','MarkerSize',11,'Linewidth',1.25);
plot(r_,d_H_levsort,styles(5),...
    'DisplayName','Levsort',...
    'MarkerSize',11,'Linewidth',1.25);
plot(r_,d_H_one_step,styles(6),...
    'DisplayName','One-step',...
    'MarkerSize',11,'Linewidth',1.25);
plot(r_,d_H_alt_min,styles(4),...
    'DisplayName','Poposed',...
    'MarkerSize',11,'Linewidth',1.25);
xticks = r_;
set(gca, 'XTick', xticks, 'XTickLabel', xticks,'Fontsize',9);
grid('on');
xlabel('$r$','interpreter','Latex','Fontsize',14);
ylabel('$d_H/n$','interpreter','Latex','Fontsize',14)
Lgnd =  legend('show');
set(Lgnd, 'Interpreter','Latex','Fontsize',12,'Location','Northwest')
title(['$n = $ ',num2str(n), ' $ m = $ ', num2str(m), ' $ d = $ ', num2str(d),...
        ' SNR $ = $' , num2str(SNR)],...
        'interpreter','Latex','Fontsize',15)
set(gca,'FontSize',18)
%ax = gca;
%exportgraphics(ax,'benchmarks.png','Resolution',300) 
%saveas(gcf,'benchmarks.fig')
