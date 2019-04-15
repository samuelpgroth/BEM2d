function[f] = basis_func_2(x_gauss,basis)
% basis is a structure

DEG = basis.deg;  % degree of polynomial
PHASE = basis.phase;

%  x1 = basis.where-basis.h/2;
%     x2 = x1+basis.h;
%     s = 0.5*(x2-x1)*x_gauss + 0.5*(x1+x2);  % map x_gauss\in[-1,1] to s\in[x1,x2]

if PHASE==0
    f = sqrt((2*DEG+1)/2).*basis.const.*legpoly(DEG+1,x_gauss);
elseif basis.which==2 % THIS SECTION IS SPECIFIC TO TRIANGLE AND NEEDS GENERALISING
    x1 = basis.where-basis.h/2;
    x2 = x1+basis.h;
    s = 0.5*(x2-x1)*x_gauss + 0.5*(x1+x2);  % map x_gauss\in[-1,1] to s\in[x1,x2]
%     f = exp(1i*PHASE*(sqrt(4*pi^2+s.^2-2*pi*s)-sqrt(3)*pi)).*sqrt((2*DEG+1)/2).*basis.const.*legpoly(DEG+1,x_gauss); 
f=sqrt((2*DEG+1)/2).*basis.const.*legpoly(DEG+1,x_gauss);
elseif basis.which==1 % the function is oscillatory
%     % Obtain s from x_gauss
    x1 = basis.where-basis.h/2;
    x2 = x1+basis.h;
    s = 0.5*(x2-x1)*x_gauss + 0.5*(x1+x2);
    f = exp(1i*PHASE*s).*sqrt((2*DEG+1)/2).*basis.const.*legpoly(DEG+1,x_gauss);
% f = (cos(PHASE*s)+sin(PHASE*s)).*sqrt((2*DEG+1)/2).*basis.const.*legpoly(DEG+1,x_gauss);
elseif basis.which==1.1
    x1 = basis.where-basis.h/2;
    x2 = x1+basis.h;
    s = 0.5*(x2-x1)*x_gauss + 0.5*(x1+x2);
    f = exp(1i*PHASE*(s-2*pi)).*sqrt((2*DEG+1)/2).*basis.const.*legpoly(DEG+1,x_gauss);
end


% keyboard

