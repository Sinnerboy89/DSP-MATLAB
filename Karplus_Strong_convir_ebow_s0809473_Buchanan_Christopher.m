clear all;

% this script further simulates an acoustic guitar string by including the
% effect of guitar body's presence (i.e. incorporating its IR). IR taken
% from a Google Drive owned by Alec Lee: https://drive.google.com/file/d/0B1XgNa5vH3j1ZkVoTFJLQ241WTg/view
% the guitar in question is a Taylor 314CE acoustic, the mic used was a
% Neumann U87. String is excited via an eBow harmonic string driver simulation, rather
% than a pluck.

Fs = 44100;
f = 110;
M = 160000; % note length
P = M/4; % eBow 'set-down' time
rho = 1;

% error checking

if rho>1
    error('rho needs to be between 0 and 1.')
elseif rho<0
    error('rho needs to be between 0 and 1.')
end

% delay line length in samples

N = round((Fs/f)-0.5);

% eBow excitation

u = KS_s0809473_Buchanan_Christopher(f,M,0.97,0.9);

x = ((1-0.99)*u); 
x(2:end) = ((1-0.99)*u(2:end))+0.99*x(1:end-1);

% hold on
% plot(x)

% Karplus-Strong-eBow algorithm

y = zeros(N+M,1);
y(1:N) = x(1:N);
y(N+1) = x(1)+rho*((y(1)/2));
for i=N+2:M+N;
    y(i) = 0.9*y(i-N)+0.1*y(i-N-1);
end

%eBow 'set-down' amplitude ramp up

rampup = linspace(0,1,P)';
y(1:P) = y(1:P).*rampup;

% output plots

% hold off
% plot(y)

% hear dry

% soundsc(y,Fs);

% convolution of dry output y with Taylor 314ce IR

ir = audioread('Taylor 314ce - Neumann U87.wav');
yc = conv(y,ir);

% hear results

soundsc(yc,Fs);