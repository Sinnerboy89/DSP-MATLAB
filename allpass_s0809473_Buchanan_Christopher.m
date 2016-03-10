function [ y ] = allpass_s0809473_Buchanan_Christopher(x,R,theta)
%ALLPASS This is a function to implement a non-time varying 2nd order all
%pass filter where x is the input signal vector, Fs is the sample rate and
%the poles are at (R,theta) and (R,-theta).

%Get recursion coefficients from pole positions:
a1 = -2*R*cos(theta);
a2 = R^2;

N = length(x);
y = zeros(N,1);

%initial conditions
y(1)= a2*x(1);
y(2)= a2*x(2)+a1*x(1)-a1*y(1);

for n=3:N
    y(n) = a2*x(n) + a1*x(n-1) + x(n-2) - a1*y(n-1) - a2*y(n-2);
end

end

