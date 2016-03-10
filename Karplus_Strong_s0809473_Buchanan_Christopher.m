clear all;

% |----------Parameters---------|

Fs = 44100;
Nyq = Fs/2;
f = 110;
M = 40000;
rho = 0.97;
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

u = (rand(N,1)*2)-1; % white noise

x = ((1-R)*u); 
x(2:end) = ((1-R)*u(2:end))+R*x(1:end-1); % filtered white noise

% |----------Time-domain plots to QC dynamics filter----------|

% plot(u)
% hold on
% plot(x)

% |----------Karplus-Strong algorithm----------|

y = zeros(N+M,1);
y(1:N) = x;
y(N+1) = x(1)+rho*((y(1)/2));
for i=N+2:M+N;
    y(i) = rho*((y(i-N)+y(i-N-1))/2);
end

% |----------Normalization----------|

y = y/(norm(y,Inf));

% |----------Output plots----------|

% hold off
% plot(y)

% l = length(y);
% fAxis = Fs*(0:l-1)'/l;
% plot(fAxis,abs(fft(y)))
% title('Frequency response')
% xlabel('frequency (Hz)')
% ylabel('|y|')
% xlim([0 Nyq])

% |----------Hear----------|

soundsc(y,Fs);

% |----------Output as 16-bit 44.1kHz wav----------|

audiowrite(['KS','.wav'],y,Fs,'BitsPerSample',16);
