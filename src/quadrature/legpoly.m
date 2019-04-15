function z=legpoly(p,x)

% input:  matrix x of points, with values \subset (-1,1)
%         p-1 the order of the Legendre polynomials (p=1,2,\ldots)
% output: a matrix containing the values of L_{p-1}(x)

%TOL=1e-12;
%if (min(min(x))+TOL<-1)|(max(max(x))-TOL>1) %added by SL 30/10/08
%    disp('warning - bad points in legpoly')
%end

z0=ones(size(x)); %L_0(x)
z1=x;             %L_1(x)
if p==1
    z=z0; %if p=1 then L_0 
elseif p==2
    z=z1; %if p=2 then L_1
else
    for j=2:p-1
        z=((2*j-1)*x.*z1-(j-1)*z0)/j; %define L_{p-1} recursively
        z0=z1;z1=z; %update L_{j-2} and L_{j-1}
    end
end
%if max(max(abs(z)))>(1e2) %ADDED 30/10/08
%    disp('max z = ')
%    disp(max(abs(z)))
%    z=(z>1)-(z<-1)+z.*(1-(z<-1)-(z>1)); %set z=-1 if z<-1 and z=1 if z>1
%end