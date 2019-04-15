function[SLout] = SLpot(x,y,w,k,q)
% function[SLout] = SLpot(x,y,w,k,q)
%
% Evaluate the single-layer potential
% x: points in domain
% y: points around the boundary
% w: quadrature weights
% k: wavenumber
% q: boundary density

x_t = size(x);
leng_x = x_t(1);
y_t = size(y);
leng_y = y_t(1);

X = repelem(x,leng_y,1);
Y = repmat(y,leng_x,1);
kernel = 1i/4*besselh(0,1,k*sqrt((X(:,1)-Y(:,1)).^2 + (X(:,2)-Y(:,2)).^2));
ko = reshape(kernel,leng_y,leng_x).';
SLout  = ko*(w.*q);