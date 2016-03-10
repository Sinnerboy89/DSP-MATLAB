%|----------PREAMBLE----------|

M = 5000; % lookback index
g = 0.9; % flange strength

%|----------AUDIO READ-IN----------|

[x, Fs]=audioread('CF_clip.wav');
if size(x,2)>1
    x = x(:,1);
end
Nyq = Fs/2;
L = length(x);
f = Fs*(0:(L/2))/L;

%|----------FEEDBACK COMB FILTER----------|

y = x;
for i=M+1:L
    y(i) = x(i) - g*y(i-M); % iteration now implicit, granting controllable feedback aspect
end

%|---------DFTs----------|

XF = fft(x);
XFM2sided = abs(XF/L);
XFM1sided = XFM2sided(1:L/2+1);
XFM1sided(2:end-1) = 2*XFM1sided(2:end-1);
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

