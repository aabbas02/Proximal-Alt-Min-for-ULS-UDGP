% download  MNIST "train set" from https://pjreddie.com/projects/mnist-in-csv/
clc
close all;
MC              = 1;
m               = 1;
n_              = [784];
r_              = [28];
options         = optimoptions('linprog','Display','none');
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
               plot_img(Y)               %title('Straight Line')
               Y_permuted      = pi_*Y;
               %title('Scrambled')
               axes(ha(2));
               plot_img(Y_permuted)               
               [energy,pi_lp]    = lp_ls_alt_min_prox(B,Y_permuted,r,lbd);
               d_H              = sum(sum(pi_ ~= pi_lp))/2;
               %title('Reconstructed')
               axes(ha(3));
               plot_img(pi_lp'*Y_permuted)
        end
    end
end
p_SNR = 10*log10(1/ ( norm(Y/max(Y) - pi_lp'*Y_permuted/max(Y_permuted),'fro')^2/n))
fig = gcf;
exportgraphics(fig,'mnist.png','Resolution',300) 
function [data,B] = get_pca_col_yale(r)
    [data] = load_yale_compressed;
    data = data(:,1:64);
    [B,~,~]=svd(data,'econ');
    B = B(:,1:r);
end
function plot_img(img)
   % imshow(reshape(round(img),[192,168]),[0,255]);
   imshow(reshape(round(img),[28,28]),[0,255]);
end
function [Y] = get_mnist()
    load('mnist_train.csv')
    Y = mnist_train(:,2:785); % exclude labels
    Y = Y';
    %load('mnist_03467.mat');
    %data = data(1:5000,:);
    %Y = data';
    %Y = Y*255;
end


function [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
% tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width 
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins 
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins 
%
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   subplot
%        pos    positions of the axes objects
%
%  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% Pekka Kumpulainen 21.5.2012   @tut.fi
% Tampere University of Technology / Automation Science and Engineering
if nargin<3; gap = .02; end
if nargin<4 || isempty(marg_h); marg_h = .05; end
if nargin<5; marg_w = .05; end

if numel(gap)==1; 
    gap = [gap gap];
end
if numel(marg_w)==1; 
    marg_w = [marg_w marg_w];
end
if numel(marg_h)==1; 
    marg_h = [marg_h marg_h];
end

axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;

py = 1-marg_h(2)-axh; 

% ha = zeros(Nh*Nw,1);
ii = 0;
for ih = 1:Nh
    px = marg_w(1);
    
    for ix = 1:Nw
        ii = ii+1;
        ha(ii) = axes('Units','normalized', ...
            'Position',[px py axw axh], ...
            'XTickLabel','', ...
            'YTickLabel','');
        px = px+axw+gap(2);
    end
    py = py-axh-gap(1);
end
if nargout > 1
    pos = get(ha,'Position');
end
ha = ha(:);
end