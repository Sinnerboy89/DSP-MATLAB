clear all;

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

% guitar chord simulation via 12 KS loops

% some example chords (input one into f).
Am = [0,110,164.8,220,261.6,329.6];
G = [98,123.5,146.8,196,246.9,392];
C = [0,130.8,164.8,196,261.6,329.6];
Em = [82.4,123.5,164.8,196,246.9,329.6];

% |----------Parameters---------|

Fs = 44100;
f = Em; % lower octave note frequency
mu = 0.2; % mu is fraction of string bridge and pluck point represents (the smaller the number, the closer to the bridge)
Ch = 0.2; % delay in seconds between each string pluck in chord (lowest to highest)

if f(1)~=0
    y1 = KS_12str_s0809473_Buchanan_Christopher(f(1),mu);
else
    y1 = 0;
end
if f(2)~=0
    y2 = KS_12str_s0809473_Buchanan_Christopher(f(2),mu);
    y2 = [zeros(round(Ch*Fs),1);y2];
else
    y2 = 0;
end
if f(3)~=0
    y3 = KS_12str_s0809473_Buchanan_Christopher(f(3),mu);
    y3 = [zeros(round(2*Ch*Fs),1);y3];
else
    y3 = 0;
end
if f(4)~=0
    y4 = KS_12str_s0809473_Buchanan_Christopher(f(4),mu);
    y4 = [zeros(round(3*Ch*Fs),1);y4];
else
    y4 = 0;
end
if f(5)~=0
    y5 = KS_12str_s0809473_Buchanan_Christopher(f(5),mu);
    y5 = [zeros(round(4*Ch*Fs),1);y5];
else
    y5 = 0;
end
if f(6)~=0
    y6 = KS_12str_s0809473_Buchanan_Christopher(f(6),mu);
    y6 = [zeros(round(5*Ch*Fs),1);y6];
else
    y6 = 0;
end

y1 = [y1;zeros(length(y6)-length(y1),1)];
y2 = [y2;zeros(length(y6)-length(y2),1)];
y3 = [y3;zeros(length(y6)-length(y3),1)];
y4 = [y4;zeros(length(y6)-length(y4),1)];
y5 = [y5;zeros(length(y6)-length(y5),1)];
y = y1+y2+y3+y4+y5+y6;

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

audiowrite(['KS_pluck_pos',mustr,'_12str_convir_chord','.wav'],y,Fs,'BitsPerSample',16);