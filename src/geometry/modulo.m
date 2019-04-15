function[remainder] = modulo(m,N)
% [remainder] = modulo(m,N)
% Different to inbuilt MATLAB function mod in that modulo(N,N) = N whereas
% mod(N,N) = 0
% S.P.Groth 30/05/12

remainder = mod(m,N);
if remainder == 0
    remainder = N;
end