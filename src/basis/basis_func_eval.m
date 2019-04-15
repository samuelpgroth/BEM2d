% basis_func_eval.m evaluates the basis function fed in at the points
% x_gauss. x_gauss are in [-1,1] and must be transformed onto [x1,x2] which
% is the support of the basis function.

function[f] = basis_func_eval(x_gauss,basis,P,tv)
% basis is a structure

DEG = basis.deg;  % degree of polynomial
PHASE = basis.phase;

x1 = basis.where-basis.h/2;
x2 = x1+basis.h;
s = 0.5*(x2-x1)*x_gauss + 0.5*(x1+x2);  % map x_gauss\in[-1,1] to s\in[x1,x2]

if PHASE==0
    f = sqrt((2*DEG+1)/2).*basis.const.*legpoly(DEG+1,x_gauss);
elseif basis.which==2 % THIS SECTION IS SPECIFIC TO TRIANGLE AND NEEDS GENERALISING
    
    f = exp(1i*PHASE*(sqrt(4*pi^2+s.^2-2*pi*s)-sqrt(3)*pi)).*sqrt((2*DEG+1)/2).*basis.const.*legpoly(DEG+1,x_gauss);

elseif basis.which==1 % the function is oscillatory
    %     % Obtain s from x_gauss
    f = exp(1i*PHASE*(s-basis.where)).*sqrt((2*DEG+1)/2).*basis.const.*legpoly(DEG+1,x_gauss);
%     f = sqrt((2*DEG+1)/2).*basis.const.*legpoly(DEG+1,x_gauss);
% f = exp(1i*PHASE*s).*sqrt((2*DEG+1)/2).*basis.const.*legpoly(DEG+1,x_gauss);

elseif basis.which==1.1
    f = exp(1i*PHASE*(s-basis.where)).*sqrt((2*DEG+1)/2).*basis.const.*legpoly(DEG+1,x_gauss);  
%  f = sqrt((2*DEG+1)/2).*basis.const.*legpoly(DEG+1,x_gauss);
% f = exp(1i*PHASE*s).*sqrt((2*DEG+1)/2).*basis.const.*legpoly(DEG+1,x_gauss);  

elseif basis.which==10
    side = basis.side;
    corner = basis.corner;
    x = repmat(P(side,:),length(s),1)+[s*tv(side,1) s*tv(side,2)];
    foo = x-repmat(P(corner,:),length(x_gauss),1);
    rtes = sqrt(foo(:,1).^2+foo(:,2).^2);
% if abs(corner-side)==2
%     r=sqrt(4*pi^2+(2*pi-s).^2);
% else
%     r=sqrt(4*pi^2+s.^2);
% end
% if norm(r-rtes)>1e-12
%     keyboard
% end
renormalise = norm(P(corner,:)-(P(side,:)+tv(side,:)*pi));
%     f = exp(1i*PHASE*(rtes-sqrt(3)*pi)).*sqrt((2*DEG+1)/2).*basis.const.*legpoly(DEG+1,x_gauss);
f = exp(1i*PHASE*(rtes-renormalise)).*sqrt((2*DEG+1)/2).*basis.const.*legpoly(DEG+1,x_gauss);

%     keyboard
%     f = legpoly(DEG+1,x_gauss);
%     f = sqrt((2*DEG+1)/2).*basis.const.*legpoly(DEG+1,x_gauss);
%     r_test = sqrt(4*pi^2+s.^2);
%     r_test2 = sqrt(4*pi^2+(2*pi-s).^2);

% % % Testing
% diffy = x-repmat(P(corner,:),length(x_gauss),1);
% dist = sqrt(diffy(:,1).^2+diffy(:,2).^2);
% % keyboard
% for i=1:length(dist)
%     direc_temp = diffy(i,:)./dist(i);
%     % Calculate refracted ray directions
%     N1 = 2; K1 = 0.025;
%     n = [tv(side,2) -tv(side,1)];
%     [N2,K2] = adjust(N1,K1,1,0,direc_temp,direc_temp,tv(side,:));
% %     keyboard
%     [d_r,direc(i,:),d_beta,d_gamma] = directions(N1,K1,N2,K2,direc_temp,direc_temp,...
%         n,tv(side,:),0,1);
% end
% 
% 
% DOT = dot(x',direc'); 
% f = exp(1i*PHASE*DOT').*sqrt((2*DEG+1)/2).*basis.const.*legpoly(DEG+1,x_gauss);
% keyboard
elseif basis.which==20
    f = exp(1i*PHASE*dot(basis.corner,tv(basis.side,:))*s).*sqrt((2*DEG+1)/2).*basis.const.*legpoly(DEG+1,x_gauss);

end

    




