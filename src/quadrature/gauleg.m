% GAULEG  computes Gauss-Legendre quadrature points 
%         suitable for integrating a smooth function 
%         on the interval [-1,1]. Can evalutate integral of a polynomial
%         of order 2n-1 with n weights and nodes exactly. Oscillations
%         will result in error. One approach is to subdivide the
%         integral into each half-period. Alternatively try methods
%         designed for oscillatory integrals, e.g. 'Filon quadrature' or 
%         'method of steepest descent'.
%
% function [x,w]=gauleg(n)
%
%input:  n = number of quadrature points
%ouput:  x = vector containing the quadrature points 
%        w = vector containing the corresponding quadrature weights
%
% Last changed by: Anne Reinarz, 2.7.2014
%
function [x,w]=gauleg(n)
  if (n < 1)
	error('number of quadrature points must be greater than 0');
  end
  x=zeros(n,1);
  w=zeros(n,1);
  m=(n+1)/2;
  xm=0.0;
  xl=1.0;
  for i=1:m
     z=cos(pi*(i-0.25)/(n+0.5));
     while 1  % run until break condition
       p1=1.0;
       p2=0.0;
       for j=1:n
          p3=p2;
          p2=p1;
          p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
       end
       pp=n*(z*p1-p2)/(z*z-1.0);
       z1=z;
       z=z1-p1/pp;
       if (abs(z-z1)<eps), break, end
     end
     x(i)=xm-xl*z;
     x(n+1-i)=xm+xl*z;
     w(i)=2.0*xl/((1.0-z*z)*pp*pp);
     w(n+1-i)=w(i);
  end
end
