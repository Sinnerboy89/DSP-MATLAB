function [ y ] = KS_12str_s0809473_Buchanan_Christopher(f,mu)
%KARPLUS_STRONG 12-STRING Frequency of note, position of pluck
%   Detailed explanation goes here

% this script further simulates an acoustic guitar string by including the
% effect of guitar body's presence (i.e. incorporating its IR), as well as picking position. IR taken
% from a Google Drive owned by Alec Lee: https://drive.google.com/file/d/0B1XgNa5vH3j1ZkVoTFJLQ241WTg/view
% the guitar in question is a Taylor 314CE acoustic, the mic used was a
% Neumann U87.

% The aim here is to simulate the sound of a 12-string, via two parallel KS
% loops, generating two different notes seperated by an octave, with
% nearly identical onset times (timing determined by D). Crude
% approximation of stiffness effect of strings is implemented by adding a
% few cents to upper octave frequency and reducing its low pass strength

Fs = 44100;
D = 0.01; % delay in seconds between first and second string being plucked
M = 160000; % initial note length
rho = 0.98;
R = 0.9; % velocity of initial pluck (~0 is a Keith Richards (NOT PETE TOWNSEND) windmill strike)
% mu is fraction of string bridge and pluck point represents (the smaller the number, the closer to the bridge)

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

if mu>1
    error('mu needs to be between 0 and 1.')
elseif mu<0
    error('mu needs to be between 0 and 1.')
end

% delay line length in samples

N1 = round((Fs/f)-0.5);
N2 = round((Fs/(2*f+round(f/45)))-0.5);

% dynamics filter

u = audioread('pluck.wav');

u(M:end) = [];

x = ((1-R)*u); 
x(2:end) = ((1-R)*u(2:end))+R*x(1:end-1);

% comb feedforward filter to simulate picking position. 

x(mu*(N1)+1:end) = x(mu*(N1)+1:end)-x(1:end-mu*(N1));

% Karplus-Strong algorithm

y1 = zeros(N1+M,1);
y1(1:N1) = x(1:N1);
y1(N1+1) = x(1)+rho*((y1(1)/2));
for i=N1+2:M+N1;
    y1(i) = rho*((y1(i-N1)+y1(i-N1-1))/2);
end
y2 = zeros(N2+M,1);
y2(1:N2) = x(1:N2);
y2(N2+1) = x(1)+rho*((y2(1)/2));
for i=N2+2:M+N2;
    y2(i) = (rho+((1-rho)/2))*((y2(i-N2)+y2(i-N2-1))/2);
end

y2 = [zeros(D*Fs,1);y2];

d = abs(length(y1)-length(y2));

if length(y1)<length(y2)
    y1 = [y1;zeros(d,1)];
elseif length(y2)<length(y1)
    y2 = [y2;zeros(d,1)];
end

y = y1+y2;

ir = audioread('Taylor 314ce - Neumann U87.wav');
y = conv(y,ir);

y = y/(norm(y,Inf));

end

