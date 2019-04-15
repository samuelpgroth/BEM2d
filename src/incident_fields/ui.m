% function[u] = ui(x,k,d,A)
% Evaluates a plane wave with
% position:   x
% wavenumber: k
% direction:  d
% amplitude:  A

function[u] = ui(x,k,d,A)

u = A*exp(1i*k*(x(:,1)*d(1)+x(:,2)*d(2)));
