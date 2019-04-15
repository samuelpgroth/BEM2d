function[SLout] = SLpotSing(xmy,k)
% function[SLout] = SLpotSing(x,y,w,k,q)
%
% Version of SLpot for singular quadrature routines
%
% Evaluate the single-layer potential
% x: points in domain
% y: points around the boundary
% w: quadrature weights
% k: wavenumber
% q: boundary density

R = abs(xmy);
SLout = 1i/4*besselh(0,1,k*R);
% keyboard
% ko = reshape(kernel,leng_y,leng_x).';
% SLout  = ko*(w.*q);