function [ y ] = feedforcomb_s0809473_Buchanan_Christopher(x,b,M)
%FEEDFORCOMB This is a function to implement an Mth order feed forward
% comb filter with coefficient b on an input x

N = length(x);
y = zeros(N,1);

%initial conditions
y(1:M) = x(1:M);

y(M+1:end) = x(M+1:end)+b*x(1:end-M);

end

