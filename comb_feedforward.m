%|----------PREAMBLE----------|

M = 5000; % lookback index

%|----------AUDIO READ-IN----------|

[x, Fs]=audioread('CF_clip.wav');
if size(x,2)>1
    x = x(:,1); % in general, input can be mono or stereo
end
L = length(x);
Nyq = Fs/2; % determines limit of certain plots
f = Fs*(0:(L/2))/L; % frequency range of input in Hz - needed for plots

%|----------FEEDFORWARD COMB FILTER----------|

y = x;
for i=M+1:L % filter can only begin to be applied once delayed signal exists
    y(i) = x(i) + x(i-M);
end

%|---------DFTs----------|

XF = fft(x);
XFM2sided = abs(XF/L); % magnitude of frequency components, positive and negative
XFM1sided = XFM2sided(1:L/2+1); % we only need see positive portion
XFM1sided(2:end-1) = 2*XFM1sided(2:end-1); % correction s.t. amplitudes are 'true' to conventional notion of positive frequency
YF = fft(y);
YFM2sided = abs(YF/L);
YFM1sided = YFM2sided(1:L/2+1);
YFM1sided(2:end-1) = 2*YFM1sided(2:end-1);


%|----------PLOTS----------|

hold on
subplot(2,2,1)
plot(f,XFM1sided)
title('Amplitude spectrum of x')
xlabel('f (Hz)')
ylabel('|x|')
xlim([0 Nyq])
subplot(2,2,3)
plot(f,XFM1sided)
title('Amplitude spectrum of x (upto 2kHz)')
xlabel('f (Hz)')
ylabel('|x|')
xlim([0 2000])
subplot(2,2,2)
plot(f,YFM1sided)
title('Amplitude spectrum of y')
xlabel('f (Hz)')
ylabel('|y|')
xlim([0 Nyq])
subplot(2,2,4)
plot(f,YFM1sided)
title('Amplitude spectrum of y (upto 2kHz)')
xlabel('f (Hz)')
ylabel('|y|')
xlim([0 2000])
