clear all;

% this script further simulates an acoustic guitar string by including the
% effect of guitar body's presence (i.e. incorporating its IR). IR taken
% from a Google Drive owned by Alec Lee: https://drive.google.com/file/d/0B1XgNa5vH3j1ZkVoTFJLQ241WTg/view
% the guitar in question is a Taylor 314CE acoustic, the mic used was a
% Neumann U87.

% |----------Parameters---------|

Fs = 44100;
f = 110;
M = 160000; % initial note length
rho = 0.98;
R = 0.9; % velocity of initial pluck (~0 is a Keith Richards (NOT PETE TOWNSEND) windmill strike)

% |----------Error checking----------|

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

% |----------Delay line length in samples----------|

N = round((Fs/f)-0.5);

% |----------Dynamics filter----------|

u = audioread('pluck.wav');

u(M:end) = [];

x = ((1-R)*u); 
x(2:end) = ((1-R)*u(2:end))+R*x(1:end-1);

% |----------Time-domain plots to QC dynamics filter----------|

% plot(u)
% hold on
% plot(x)

% |----------Karplus-Strong algorithm----------|

y = zeros(N+M,1);
y(1:N) = x(1:N);
y(N+1) = x(1)+rho*((y(1)/2));
for i=N+2:M+N;
    y(i) = rho*((y(i-N)+y(i-N-1))/2);
end

% |----------Output plots----------|

% hold off
% plot(y)

% |----------Convolution of dry output y with Taylor 314ce IR----------|

ir = audioread('Taylor 314ce - Neumann U87.wav');
yc = conv(y,ir);

% |----------Normalization---------|

y = y/(norm(y,Inf));
yc = yc/(norm(yc,Inf));

% |----------Compare dry and wet in time domain----------|

% hold on
% plot(yc)
% plot(y)

% |----------Hear results----------|

soundsc(y,Fs);
pause(2);
soundsc(yc,Fs);

% |----------Output as 16-bit 44.1kHz wav----------|

audiowrite(['KS_pluck_convir','.wav'],yc,Fs,'BitsPerSample',16);