close all
clc
clear all
m    = 1;
SNR_ = [0:10:100];
table = zeros(2*2,1*length(SNR_));
for i = 1 : length(SNR_)
    SNR = SNR_(i);
    n   = 5000;
    r  = n/20;
    k  = n/2;
    H_r = stirling_function(r)*(n/r);
    H_k = stirling_function(n) - stirling_function(k);
    table(4,i) = 1 - ( 1 + (m*n/2)*(log2(1+SNR) ) )/H_r;
    table(3,i) = 1 - ( 1 + (m*n/2)*(log2(1+SNR) ) )/H_k;
    n   = 1000;
    r  = n/10;
    k  = n/2;
    H_r = stirling_function(r)*(n/r);
    H_k = stirling_function(n) - stirling_function(k);
    table(2,i) = 1 - ( 1 + (m*n/2)*(log2(1+SNR) ) )/H_r;
    table(1,i) = 1 - ( 1 + (m*n/2)*(log2(1+SNR) ) )/H_k;
end
figure;
hold on

plot(SNR_,table(2,:),'r-*',...
    'DisplayName','r-local, $r = n/10$, $n=1000$',...
     'MarkerSize',9,'Linewidth',1.75);
plot(SNR_,table(1,:),'b-*',...
    'DisplayName','k-sparse, $k = n/2$',...
     'MarkerSize',9,'Linewidth',1.75); 
 plot(SNR_,table(4,:),'r--x',...
    'DisplayName','r-local, $r = n/20$, $n=5000$',...
     'MarkerSize',9,'Linewidth',1.75);
plot(SNR_,table(3,:),'b--x',...
    'DisplayName','k-sparse, $k = n/2$',...
     'MarkerSize',9,'Linewidth',1.75); 
set(gca, 'XTick', xticks, 'XTickLabel', xticks,'Fontsize',14);
set(gca, 'YTick', yticks, 'YTickLabel', yticks,'Fontsize',14);
grid('on');
ylabel('$\Pr\{{\hat{\mathbf{P}} \neq \mathbf{P}}\}$','interpreter','Latex','Fontsize',16);
xlabel('necessary SNR','interpreter','Latex','Fontsize',16)
Lgnd =  legend('show');
set(Lgnd, 'Interpreter','latex','Location','northeast','Fontsize',16)
ax = gca;
%exportgraphics(ax,'nec_snr.png','Resolution',300) 
% https://en.wikipedia.org/wiki/Stirling%27s_approximation - Speed of convergence and error estimates
function x =  stirling_function(n)
x = n*log(n) - n + 0.5*log(2*pi*n) + (1/(12*n)) - (1/(360*n^3)) + 1/(1260*n^5) - 1/(1680*n^7);
x = x/log(2);
end