function [ x,z,r, w] = DuffyLog( a, b, N, width_A)
%Designed for functions on S=[a,b]x[a,b] with logarithmic singularity on diagonal {(x,y)\in S:
%x=y}. Returns N^2 generalised gauss points in S. Returns the distance r to
%avoid rounding errors.
  
    if nargin==3
        width_A=b-a; %let the width_A input be optional.
    end
    
    %compute weights and nodes for L[f]
    [x,w] = quad_gengauss_log(N);
    [gamma,mu] = rearrange(x,x);
    [w1,w2] = rearrange(w,w);
    w_L=w1.*w2;
    %weights and nodes over R[f] will start identically
    w_R=w_L;	sigma=gamma;   zeta=mu;

    W_L=width_A^2.*mu.*w_L;
    x_L=a+mu*width_A;
    z_L=a+(1-gamma).*mu*width_A;
    r_L=gamma.*mu*width_A;

    W_R=w_R.*zeta*width_A^2;
    x_R=a+(1-zeta)*width_A;
    z_R=width_A*(sigma.*zeta+1-zeta)+a;
    r_R=width_A*zeta.*sigma;

    w=[W_L.' W_R.'].';
    r=[r_L.' r_R.'].';
    x=[x_L.' x_R.'].';
    z=[z_L.' z_R.'].';
end