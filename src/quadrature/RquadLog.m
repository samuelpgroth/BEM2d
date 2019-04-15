function [x,y,r_,w] = RquadLog( a, b, c,  N, angle, A_width, C_width)
%This code outputs quadrature points x,y,r:=|x-y| and weights w.
%It is designed for logarithmic singularities on the corner of a square [a,b]x[b,c]
%This is acheived by mapping the square to a quarter circle, and using
%polar coordinates - with the singularity only at r=0

%the "angle" input corresponds to the angle in 2 dimensions at point b, for
%example if [a,b] parameterises one side of a triangle and [b,c]
%parameterises another. This can be ignored and will be set to pi
%(corresponding to a straight line) by default. Similarly, exact widths can
%be entered in the final two args, alternatively these will be computed
%to machine precision.

if nargin<=5
    A_width=b-a; C_width=c-b;
end

if nargin==4
    angle=pi;
end

    [r,wr] = quad_gengauss_log(N);
    [theta1,w1] = gauss_quad(pi/4,pi/2,N);
    [theta2,w2] = gauss_quad(0,pi/4,N);

    %Will need 2N repetitions of vector r to match with thetas (+ weights). But
    %actually make N repetitions, and use twice, once with each theta.
    R=repmat(r,N,1);
    Wr=repmat(wr,N,1);

    %Now create a stretched vector of each theta, with N repetitions of each
    %entry. And do the same for their weights.
    Th1=reshape(repmat(theta1.',N,1),N*N,1);
    Wt1=reshape(repmat(w1.',N,1),N*N,1);

    Th2=reshape(repmat(theta2.',N,1),N*N,1);
    Wt2=reshape(repmat(w2.',N,1),N*N,1);

    %upper quarter circle first
    xU=b-A_width*R.*cot(Th1);   yU=b+C_width*R;    
    %lower quarter circle
    xL=b-A_width*R;   yL=b+C_width*R.*tan(Th2);          
    
    %And the weights
    wU=Wt1.*Wr.*R.*csc(Th1).^2*A_width*C_width;
    wL=Wt2.*Wr.*R.*sec(Th2).^2*A_width*C_width;
    %Now stick vectors together
    x=[xU.' xL.'].';  y=[yU.' yL.'].';  w=[wU.' wL.'].';

    %also output the r_=|x-y| values

    if angle==pi
        r_=[R.*(A_width*cot(Th1)+C_width)
            R.*(A_width+C_width*tan(Th2))];
    else
        r_=[R.*((A_width*cot(Th1)).^2 + C_width.^2 - 2*cot(Th1)*A_width*C_width*cos(angle)).^(1/2)
            R.*(A_width.^2 + (C_width*tan(Th2)).^2 - 2*tan(Th2)*A_width*C_width*cos(angle)).^(1/2)];
    end

end
%-------------------------------------------------------------------------%
%Testing im maple against:

%maple:
%int(int(HankelH1(0, abs(x-y))*x^2*y^3, x = 0 .. 1), y = 1 .. 2.)
%    =.98470493593248650184-0.69944511534011848300e-1i

%matlab:
%   f=@(x,y,r)besselh(0,1,r).*x.^2.*y.^3;
%   [x,y,r,w] = RquadLog( 0, 1, 2, N, pi, 1,1);
%   f(x,y,r).'*w;
%   error=abs(f(x,y,r).'*w - (.98470493593248650184-0.69944511534011848300e-1i))

%error      N
%O(e-13)    10
%O(e-16)	50
%O(e-16)	100