clear all;

% this script further simulates an acoustic guitar string by including the
% effect of guitar body's presence (i.e. incorporating its IR), as well as picking position. IR taken
% from a Google Drive owned by Alec Lee: https://drive.google.com/file/d/0B1XgNa5vH3j1ZkVoTFJLQ241WTg/view
% the guitar in question is a Taylor 314CE acoustic, the mic used was a
% Neumann U87.

% The aim here is to simulate the sound of a mandolin, via two parallel KS
% loops, generating two different notes seperated by a few cents (determined by fd), with
% nearly identical onset times (timing determined by D).

% |----------Parameters---------|

Fs = 44100;
cf = 220; % central note frequency
fd = round(cf/220); % 2*fd = difference in frequency between two strings. 220 selected as that's when 'beating' starts to become audible
D = 0.01; % delay in seconds between first and second string being plucked
M = 160000; % initial note length
rho = 0.98;
R = 0.9; % velocity of initial pluck (~0 is a Keith Richards (NOT PETE TOWNSEND) windmill strike)
mu = 0.1; % mu is fraction of string bridge and pluck point represents (the smaller the number, the closer to the bridge)

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

N1 = round((Fs/(cf-fd))-0.5);
N2 = round((Fs/(cf+fd))-0.5);

% dynamics filter

u = audioread('pluck.wav');

u(M:end) = [];

x = ((1-R)*u); 
x(2:end) = ((1-R)*u(2:end))+R*x(1:end-1);

% comb feedforward filter to simulate picking position. 

x(mu*(round((N1+N2)/2))+1:end) = x(mu*(round((N1+N2)/2))+1:end)-x(1:end-mu*(round((N1+N2)/2)));

% time-domain plots to QC dynamics filter

% plot(u)
% hold on
% plot(x)

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
    y2(i) = rho*((y2(i-N2)+y2(i-N2-1))/2);
end

y2 = [zeros(D*Fs,1);y2];

d = abs(length(y1)-length(y2));

if length(y1)<length(y2)
    y1 = [y1;zeros(d,1)];
elseif length(y2)<length(y1)
    y2 = [y2;zeros(d,1)];
end

y = y1+y2;

% output plots

% hold off
% plot(y)

% hear dry

% soundsc(y,Fs);

% convolution of dry output y with Taylor 314ce IR

ir = audioread('Taylor 314ce - Neumann U87.wav');
yc = conv(y,ir);

yc = yc/(norm(yc,Inf));

% compare dry and wet in time domain

% hold on
% plot(yc)
% plot(y)

% hear results

soundsc(yc,Fs);

% |---------Convert useful parameters to strings for filename----------|

mustr = num2str(mu);

% |----------Output as 16-bit 44.1kHz wav----------|

audiowrite(['KS_pluck_pos',mustr,'_mand_convir','.wav'],y,Fs,'BitsPerSample',16);

