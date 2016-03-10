clear all;

% this script uses the provided pluck.wav instead of filtered white noise

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

u = audioread('pluck.wav'); % pluck sample taken from Learn

u(M:end) = []; % clip sample temporally

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

% |----------Normalization----------|

y = y/(norm(y,Inf));

% |----------Output plots----------|

% hold off
% plot(y)

% |----------Hear----------|

soundsc(y,Fs);

% |----------Output as 16-bit 44.1kHz wav----------|

audiowrite(['KS_pluck','.wav'],y,Fs,'BitsPerSample',16);