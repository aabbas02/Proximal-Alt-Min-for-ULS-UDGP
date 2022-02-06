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
d               = 200;
r               = 200;
m_              = [1 5 10 15 20 25 30 35 40 45 50];
n_              = [1000 1200 1400 1600 1800 2000 2200 2400 2600 2800 3000];
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
            X_hat                = B\(pi_lp'*Y_permuted_noisy);
            err_X(i,j)           = norm(X_hat - X,'fro')/norm(X,'fro')...
                                   + err_X(i,j);
        end
    end
end
err_X  = err_X/MC;
err_Pi = err_Pi/MC;















close all
imagesc(flip(err_X))
colormap(gray)
set(gca,'ytick',1:length(n_))
set(gca,'xtick',1:length(n_))
set(gca,'yticklabels',flip(n_),'Fontsize',12)
set(gca,'xticklabels',m_)
xlabel('$m$','Interpreter','Latex','Fontsize',15)
ylabel('$n$','Interpreter','Latex','Fontsize',15)
colorbar
%title('$||\mathbf{x} - \hat \mathbf{x}||/||\mathbf{x}||$','Interpreter','Latex','Fontsize',15) 
title(['Error in $\widehat \mathbf X$. $d = $',num2str(d),...
        ', $r = $ ',num2str(r),...
        ', $\mathrm{SNR} = $', num2str(SNR),'.'],...
        'Interpreter','Latex','Fontsize',13) 
c = colorbar;
c.Label.String = '$\frac{||\mathbf{x}^* - \hat \mathbf{x}||}{||\mathbf{x}^*||}$';
c.Label.Interpreter = 'Latex';
c.Label.Position = [pos(1)/2 pos(2)+1.20]; % to change its position
c.Label.Rotation = 0; % to rotate the text
set(gca,'FontSize',13)
ax = gca;
c.Label.FontSize = 14;
exportgraphics(ax,'rel_err_x.pdf','Resolution',300) 
saveas(gcf,'rel_err_x.fig')
figure
imagesc(flip(err_Pi))
colormap(gray)
set(gca,'ytick',1:length(n_))
set(gca,'xtick',1:length(n_))
set(gca,'yticklabels',flip(n_),'Fontsize',12)
set(gca,'xticklabels',m_)
xlabel('$m$','Interpreter','Latex','Fontsize',14)
ylabel('$n$','Interpreter','Latex','Fontsize',14)
c = colorbar;
c.Label.String = '$d_H/n$';
c.Label.Interpreter = 'Latex';
c.Label.Position = [pos(1)/2 pos(2)+0.95]; % to change its position
c.Label.Rotation = 0; % to rotate the text
title(['Error in $\widehat \mathbf P$. $d = $',num2str(d),...
        ', $r = $ ',num2str(r),...
        ', $\mathrm{SNR} = $ ', num2str(SNR),'.'],...
        'Interpreter','Latex','Fontsize',13) 
ax = gca;
set(gca,'FontSize',13)
c.Label.FontSize = 14;
exportgraphics(ax,'rel_err_P.pdf','Resolution',300) 
saveas(gcf,'rel_err_P.fig')

