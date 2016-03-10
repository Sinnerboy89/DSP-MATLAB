clear all;

[x,Fs] = audioread('Dance_Stereo.wav');
x = (x(:,1)+x(:,2))/2;
Nyq = Fs/2;

theta = 2;
R = 0.8;

y = allpass(x,R,theta);

N = length(x);
f = Fs*(0:(N/2))/N;
fAxis = Fs*(0:N-1)/N;

XF = fft(x);
YF = fft(y);
HFM = abs(YF)./abs(XF);
HFP = (angle(YF)-angle(XF));

subplot(2,1,1)
plot(fAxis,HFM)
title('Frequency response')
xlabel('frequency (Hz)')
ylabel('|H|')
xlim([0 Nyq])
subplot(2,1,2)
plot(fAxis,HFP)
title('Phase response')
xlabel('frequency (Hz)')
ylabel('ang(H)')
xlim([0 Nyq])

% figure
% subplot(2,1,1)
% plot(fAxis,HFM)
% title('Frequency response (upto 20kHz)')
% xlabel('frequency (Hz)')
% ylabel('|H|')
% xlim([0 20000])
% subplot(2,1,2)
% plot(fAxis,HFP)
% title('Phase response (upto 20kHz)')
% xlabel('frequency (Hz)')
% ylabel('ang(H)')
% xlim([0 20000])