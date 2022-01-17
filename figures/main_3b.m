str = pwd;
k   = strfind(str,'\');
str = str(1:k(end));
cd (str)
addpath('code\misc','code\prox_alt_min') 
clc
close all;
clear all
MC              = 5;
SNR             = 100;
d               = 100;
m_              = [25 50];
r_              = [10 25 50 100];
n               = 1000;
d_H_lp_col      = zeros(2,length(r_));
options         = optimoptions('linprog','Display','none');
lbd             = 100;
for j = 1 : length(r_) 
	r           = r_(j);
    for k = 1 : MC
        B           = randn(n,d);
        X           = randn(d,max(m_));
        Y_full      = B*X;  
        pi_         = make_r_local_permutation(n,r);
            for i = 1 : 2
                m                = m_(i);
                Y                = Y_full(:,1:m);
                noise_var   	 = norm(X,'fro')^2  / (SNR*m);
                W_big            = sqrt(noise_var)*randn(n,m);
                Y_permuted       = pi_*Y;
                Y_permuted_noisy = Y_permuted + W_big;
                [~,pi_lp]        = lp_ls_alt_min_prox(B,Y_permuted_noisy,r,lbd);
                d_H              = sum(sum(pi_ ~= pi_lp))/(2*n);
                d_H_lp_col(i,j)  = d_H + d_H_lp_col(i,j);
            end
    end
end
d_H_lp_col = d_H_lp_col/MC;
close all
clc
%styles =["g-x","b-x","r-x","k-x"];
hold on;
plot(r_,d_H_lp_col(1,:),"k--x",...
    'DisplayName',['$m = $',num2str(m_(1))],...
    'MarkerSize',15,'Linewidth',1.75);
plot(r_,d_H_lp_col(2,:),"k-x",...
    'DisplayName',['$m = $',num2str(m_(2))],...
    'MarkerSize',15,'Linewidth',1.75);
xticks = r_;
yticks =  0.05*ceil(d_H_lp_col/0.05);
yticks = unique(yticks);
yticks(yticks==0)=[];
set(gca, 'XTick', xticks, 'XTickLabel', xticks,'Fontsize',9);
set(gca, 'YTick', yticks, 'YTickLabel', yticks,'Fontsize',9);
grid('on');
xlabel('$r$','interpreter','Latex','Fontsize',15);
ylabel('$d_H/n$','interpreter','Latex','Fontsize',15)
Lgnd =  legend('show');
set(Lgnd, 'Interpreter','Latex','Fontsize',17,'Location','Northwest');
title(['$n = $ ',num2str(n),  ' $ d = $ ', num2str(d),...
        ' SNR $ = $' , num2str(SNR)],...
        'interpreter','Latex','Fontsize',15)
set(gca,'FontSize',18)
%ax = gca;
%exportgraphics(ax,'n_1000_d_100_r_100.png','Resolution',300) 
%saveas(gcf,'n_1000_d_100.fig')
