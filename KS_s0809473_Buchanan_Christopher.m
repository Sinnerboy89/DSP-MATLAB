function [ y ] = KS_s0809473_Buchanan_Christopher(f,M,rho,R)
% KS Output of simple KS algorithm
%   Detailed explanation goes here

Fs = 44100;

% error checking

if rho>1
    error('rho needs to be between 0 and 1.')
elseif rho<0
    error('rho needs to be between 0 and 1.')
end

if R>1
    error('R needs to be between 0 and 1.')
elseif R<0
    error('R needs to be between 0 and 1.')
end

% delay line length in samples

N = round((Fs/f)-0.5);

% dynamics filter

u = (rand(N,1)*2)-1;

x = ((1-R)*u); 
x(2:end) = ((1-R)*u(2:end))+R*x(1:end-1);


% Karplus-Strong algorithm

y = zeros(N+M,1);
y(1:N) = x;
y(N+1) = x(1)+rho*((y(1)/2));
for i=N+2:M+N;
    y(i) = rho*((y(i-N)+y(i-N-1))/2);
end

end

