function [y] = allpassTV_s0809473_Buchanan_Christopher(x,R,theta, delta,rate,Fs)
%ALLPASS This is a function to implement a non-time varying 2nd order all
%pass filter where x is the input signal vector, Fs is the sample rate and
%the poles are at (R,theta) and (R,-theta).

N = length(x);
m = 1:N;
thetavec = theta+(delta*(1-cos(2*pi*m*rate/Fs)));
y = zeros(N,1);

%Get recursion coefficients from pole positions:
a1 = -2*R*cos(thetavec);
a2 = R^2;

%initial conditions
y(1)= a2*x(1);
y(2)= a2*x(2)+a1(2)*x(1)-a1(2)*y(1);

for n=3:N
    y(n) = a2*x(n) + a1(n)*x(n-1) + x(n-2) - a1(n)*y(n-1) - a2*y(n-2);
end

end

