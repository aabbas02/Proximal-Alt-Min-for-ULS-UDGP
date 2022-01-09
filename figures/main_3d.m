str = pwd;
k   = strfind(str,'\');
str = str(1:k(end));
cd (str)
addpath('misc','prox_alt_min','udgp') 
clc
close all;
MC              = 25;
sigma_          = [-5 -4.5 -4 -3.5 -3 -2.5 -2 -1.5 -1.0];
sigma_          = 10.^(sigma_);
d               = 25;
n               = d;
m               = 1;
options         = optimoptions('linprog','Display','none');
B = zeros(n*(n-1)/2,n-1);
B(1:n-1,1:n-1) = eye(n-1);
pos1 = 1;
pos2 = 2;
lbd = 1;
for row = n : n*(n-1)/2
    B(row,pos1) = 1;
    B(row,pos2) = -1;
    if(pos2 == n-1)
        pos1 = pos1+1;
        pos2 = pos1+1;
    else
        pos2 = pos2+1;
    end
end
err = zeros(3,length(sigma_));
pi_ = get_r_local_perm_udgp(n);
for t = 1 : 3
    if(t==1)
        x           = randn(d,1);
    elseif (t==2)
        x           = sqrt(5)*randn(d,1);
    else 
        x           = sqrt(10)*randn(d,1);
    end
    x               = x - min(x);
    x               = sort(x);
    x(2:end)        = x(end:-1:2);
    Y               = B*x(2:end);
for i = 1 : length(sigma_)
    noise_var   	= sigma_(i);
    for k = 1 : MC
                W_big            = sqrt(noise_var)*randn(n*(n-1)/2,m);
                Y_permuted       = pi_*Y;
                Y_permuted_noisy = Y_permuted + W_big;
                [pi_lp]          = lp_ls_alt_min_prox_udgp(B,Y_permuted_noisy,n,lbd);
                x_hat = B\(pi_lp'*Y_permuted_noisy);
                x_hat = [0;x_hat];
                x_hat = sort(x_hat);
                x_hat(2:end) = x_hat(end:-1:2);
                err(t,i)  = err(t,i) + norm(x_hat-x)/norm(x);
                k
    end
end
end
close all
clc
styles =["g-x","b-x","r-x","k-x"];
hold on;
plot(log10(sigma_),err(1,:),styles(1),...
    'DisplayName',['$\mathbf{N}(0,1)$'],...
    'MarkerSize',15,'Linewidth',1.75);
plot(log10(sigma_),err(2,:),styles(3),...
    'DisplayName',['$\mathbf{N}(0,5)$'],...
    'MarkerSize',15,'Linewidth',1.75);
plot(log10(sigma_),err(3,:),styles(2),...
    'DisplayName',['$\mathbf{N}(0,10)$'],...
    'MarkerSize',15,'Linewidth',1.75);
xticks = log10(sigma_);
yticks = cat(2,err(1,:),err(2,:),err(3,:));
yticks =  0.2*round(yticks/0.2);
yticks = unique(yticks);
set(gca, 'XTick', xticks, 'XTickLabel', xticks,'Fontsize',9);
set(gca, 'YTick', yticks, 'YTickLabel', yticks,'Fontsize',9);
grid('on');
Lgnd =  legend('show');
set(Lgnd, 'Interpreter','Latex','Fontsize',15,'Location','Northwest');
title('Reconstruction error uDGP, $d = 50$ points','interpreter','Latex')
ylabel('$\frac{||\mathbf{x^*} - \hat{\mathbf{x}}||}{||\mathbf{x}||}$','interpreter','Latex')
set(gca,'FontSize',13)
xlabel('Noise variance $\log_{10} \sigma^2$','interpreter','Latex','Fontsize',17);
%ylabel('$\frac{||\mathbf{x^}* - \hat{\mathbf{x}}||}{||\mathbf{x}||}$','interpreter','Latex','Fontsize',19)
ax = gca;
%exportgraphics(ax,'udgp.png','Resolution',300) 
saveas(gcf,'udgp_2.fig')


