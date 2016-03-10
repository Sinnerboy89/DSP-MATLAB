clear all

[x,Fs] = audioread('DryGuitar_mono_44100_481344.wav');

if size(x,2)>1
    x = (x(:,1)+x(:,2))/2;
    disp('input is stereo - converted to mono')
end

Nyq = Fs/2;

N = length(x);
M = 10;
g = 0.5;

y = allpassM_s0809473_Buchanan_Christopher(x,g,M);

b=zeros(1,M+1);
b(1) = g;
b(M+1) = 1;

a=zeros(1,M+1);
a(1) = 1;
a(M+1) = g;

fAxis = Fs*(0:N-1)'/N;
tAxis = (0:N-1)'/Fs;

XF = fft(x);
YF = fft(y);
HF = YF./XF;
HFM = abs(YF)./abs(XF);
HFP = (angle(YF)-angle(XF));
h = ifft(HF);

subplot(2,3,1)
plot(abs(HF))
title('Frequency response')
xlabel('frequency (Hz)')
ylabel('|H|')
xlim([0 Nyq])
subplot(2,3,2)
plot(fAxis,unwrap(HFP))
title('Phase response')
xlabel('frequency (Hz)')
ylabel('ang(H)')
xlim([0 Nyq])
subplot(2,3,3)
zplane(b,a)
title('Pole/Zero plot')
subplot(2,3,4)
plot(x)
hold on
plot (y)
title('Time domain test signal')
xlabel('Time (samples)')
subplot(2,3,5)
plot(tAxis,h)
title('Impulse response h')
xlabel('Time (s)')
xlim([0 0.05])