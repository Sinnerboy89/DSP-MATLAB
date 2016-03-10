clear all;

% this script allows alpha to be frequency band dependent

Fs = 44100;
Ts = 1/Fs;
Nyq = Fs/2;

% number of reflections allowed in one dimension

N = 10;

% room dimensions

Lx = 10; 
Ly = 10;
Lz = 10;

% wall absorption coefficients

alpha1 = 0.9;
alpha2 = 0.9;
alpha3 = 0.9;
alpha4 = 0.9;
alpha5 = 0.9;
alpha6 = 0.9;
alpha7 = 0.9;
alpha8 = 0.9;
alpha9 = 0.9;

% position of source

p = 2;
q = 7;
r = 1;

% position of receiver

a = 2;
b = 4;
c = 5;

% IR calculation initializations

DELN = zeros((2*N+1),(2*N+1),(2*N+1));
G1 = DELN;
Nspan = (-N:N);
Aodd = Nspan*Lx-p-a;
Aeven = Nspan*Lx+p-a;
Bodd = Nspan*Ly-q-b;
Beven = Nspan*Ly+q-b;
Codd = Nspan*Lz-r-c;
Ceven = Nspan*Lz+r-c;
delnf = sqrt((((N+1)*Lx)^2)+(((N+1)*Ly)^2)+(((N+1)*Lz)^2))*Fs/343;
h1 = zeros(ceil(delnf),1);
h2 = h1;
h3 = h1;
h4 = h1;
h5 = h1;
h6 = h1;
h7 = h1;
h8 = h1;
h9 = h1;
ti = Fs/343;

% DELN (delay times in samples) and G (amplitudes) computation

for d=-N:N;
    for e=-N:N;
        for f=-N:N;
            if mod(d,2)==1
                A = Aodd(d+N+2);
            else
                A = Aeven(d+N+1);
            end
            if mod(e,2)==1
                B = Bodd(e+N+2);
            else
                B = Beven(e+N+1);
            end
            if mod(f,2)==1
                C = Codd(f+N+2);
            else
                C = Ceven(f+N+1);
            end
            l = sqrt((A^2)+(B^2)+(C^2));
            t = ti*l;
            w = abs(d)+abs(e)+abs(f);
            DELN(d+N+1,e+N+1,f+N+1) = round(t);
            G1(d+N+1,e+N+1,f+N+1) = ((alpha1)^w)/l;
            h1(DELN(d+N+1,e+N+1,f+N+1)) = h1(DELN(d+N+1,e+N+1,f+N+1)) + G1(d+N+1,e+N+1,f+N+1);
            G2(d+N+1,e+N+1,f+N+1) = ((alpha2)^w)/l;
            h2(DELN(d+N+1,e+N+1,f+N+1)) = h2(DELN(d+N+1,e+N+1,f+N+1)) + G2(d+N+1,e+N+1,f+N+1);
            G3(d+N+1,e+N+1,f+N+1) = ((alpha3)^w)/l;
            h3(DELN(d+N+1,e+N+1,f+N+1)) = h3(DELN(d+N+1,e+N+1,f+N+1)) + G3(d+N+1,e+N+1,f+N+1);
            G4(d+N+1,e+N+1,f+N+1) = ((alpha4)^w)/l;
            h4(DELN(d+N+1,e+N+1,f+N+1)) = h4(DELN(d+N+1,e+N+1,f+N+1)) + G4(d+N+1,e+N+1,f+N+1);
            G5(d+N+1,e+N+1,f+N+1) = ((alpha5)^w)/l;
            h5(DELN(d+N+1,e+N+1,f+N+1)) = h5(DELN(d+N+1,e+N+1,f+N+1)) + G5(d+N+1,e+N+1,f+N+1);
            G6(d+N+1,e+N+1,f+N+1) = ((alpha6)^w)/l;
            h6(DELN(d+N+1,e+N+1,f+N+1)) = h6(DELN(d+N+1,e+N+1,f+N+1)) + G6(d+N+1,e+N+1,f+N+1);
            G7(d+N+1,e+N+1,f+N+1) = ((alpha7)^w)/l;
            h7(DELN(d+N+1,e+N+1,f+N+1)) = h7(DELN(d+N+1,e+N+1,f+N+1)) + G7(d+N+1,e+N+1,f+N+1);
            G8(d+N+1,e+N+1,f+N+1) = ((alpha8)^w)/l;
            h8(DELN(d+N+1,e+N+1,f+N+1)) = h8(DELN(d+N+1,e+N+1,f+N+1)) + G8(d+N+1,e+N+1,f+N+1);
            G9(d+N+1,e+N+1,f+N+1) = ((alpha9)^w)/l;
            h9(DELN(d+N+1,e+N+1,f+N+1)) = h9(DELN(d+N+1,e+N+1,f+N+1)) + G9(d+N+1,e+N+1,f+N+1);
        end
    end
end

L = length(h1);
fh = Fs*(0:(L/2))/L;

H1F = fft(h1);
H1F = H1F(1:ceil(L/2));
keepInd = fh <= 62; % Keep frequencies between 0 and 62 Hz
H1F(~keepInd) = 0;

H2F = fft(h2);
H2F = H2F(1:ceil(L/2));
keepInd = fh >= 63 & fh <= 125; % Keep frequencies between 63 and 125 Hz
H2F(~keepInd) = 0;

H3F = fft(h3);
H3F = H3F(1:ceil(L/2));
keepInd = fh >= 126 & fh <= 250; % Keep frequencies between 126 and 250 Hz
H3F(~keepInd) = 0;

H4F = fft(h4);
H4F = H4F(1:ceil(L/2));
keepInd = fh >= 251 & fh <= 500; % Keep frequencies between 251 and 500 Hz
H4F(~keepInd) = 0;

H5F = fft(h5);
H5F = H5F(1:ceil(L/2));
keepInd = fh >= 501 & fh <= 1000; % Keep frequencies between 501 and 1000 Hz
H5F(~keepInd) = 0;

H6F = fft(h6);
H6F = H6F(1:ceil(L/2));
keepInd = fh >= 1001 & fh <= 2000; % Keep frequencies between 1001 and 2000 Hz
H6F(~keepInd) = 0;

H7F = fft(h7);
H7F = H7F(1:ceil(L/2));
keepInd = fh >= 2001 & fh <= 4000; % Keep frequencies between 2001 and 4000 Hz
H7F(~keepInd) = 0;

H8F = fft(h8);
H8F = H8F(1:ceil(L/2));
keepInd = fh >= 4001 & fh <= 8000; % Keep frequencies between 4001 and 8000 Hz
H8F(~keepInd) = 0;

H9F = fft(h9);
H9F = H9F(1:ceil(L/2));
keepInd = fh >= 8001 & fh <= Nyq; % Keep frequencies between 8001 and Nyq Hz
H9F(~keepInd) = 0;

HF = H1F+H2F+H3F+H4F+H5F+H6F+H7F+H8F+H9F;
h = ifft(HF,'symmetric');

% All reflections included check

% if delnf<max(max(max(DELN)))
%     disp ('DELN bad...somewhere')
% elseif delnf>max(max(DELN))
%     disp ('DELN all good.')
% end

% plot DELN as a spatial QC. Animation moves through
% dimension varied by i

% for i=1:21
%     stem3(-N:N,-N:N,squeeze(DELN(i,:,:)))
%     drawnow
%     pause(0.1)
% end

% normalize IR

h = h/(norm(h,Inf));
l = length(h);

% plot impulse response

T = (0:Ts:(l-1)*Ts)';
plot(T,h)
xlabel('time (s)')
ylabel('amplitude')
title('IR')

% plot frequency response

% fAxis = Fs*(0:l-1)'/l;
% HFM = abs(fft(h));
% plot(HFM)
% title('Frequency response')
% xlabel('frequency (Hz)')
% ylabel('|H|')
% xlim([0 Nyq])

%generate impulse response filename

Lxstr = num2str(Lx);
Lystr = num2str(Ly);
Lzstr = num2str(Lz);
alphastr = num2str(alpha1);
Nstr = num2str(N);
params = strcat(Lxstr,'x',Lystr,'x',Lzstr,'_','_freqalpha','_N',Nstr);

%output IR as 16-bit 44.1kHz wav

audiowrite(['ir_',params,'.wav'],h,Fs,'BitsPerSample',16);




