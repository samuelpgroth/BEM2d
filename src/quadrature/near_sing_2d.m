function[x,y,r,w] = near_sing_2d(a,b,c,d,xG,wG,n_layers)
% [x,y,r,w] = near_sing_2d(a,b,c,d,xG,wG,n_layers)
%
% Sets up graded quadrature for the near singularity located. Only call
% this if |c-b|<0.15, i.e. there is a singularity near domain of
% integration.
% 
% [a,b]    : domain of integration in x
% [c,d]    : domain of integration in y
% xG,wG    : standard Gaussian nodes in [-1,1] and weights
% n_layers : number of layers in grading

if a>d  % then switch [a,b] and [c,d]
    C = a;      D = b;
    A = c;      B = d;
else
    A=a;B=b;C=c;D=d;
end

s = (B*D-A*C)/(D-C+B-A);   % "geometric middle" where singularity is

% Shift a,c,b,d so that we are grading towards 0 for round-off error
% purposes.
A = A-s; B = B-s;
C = C-s; D = D-s;

[x_grade,w_grade] = near_sing_1d(-B,-A,-B,xG,wG,n_layers);
x_grade = -flipud(x_grade); % since we the grading is always done to the left in near_sing_1d
w_grade = flipud(w_grade);
[y_near,w_near] = near_sing_1d(C,D,C,xG,wG,n_layers);

% Shift back
x_grade = x_grade+s;
y_near = y_near+s;

[z,w] = TensorQuad(x_grade,w_grade,y_near,w_near);


if a>d  % then switch [a,b] and [c,d]
    x = z(:,2); y = z(:,1);
else
    x = z(:,1); y = z(:,2);
end

r = abs(x-y);




