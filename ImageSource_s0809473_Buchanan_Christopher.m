clear all;

Fs = 44100;
Ts = 1/Fs;
Nyq = Fs/2;

% number of reflections allowed in one dimension

N = 10;

% room dimensions

Lx = 10; 
Ly = 10;
Lz = 10;

% wall absorption coefficient

alpha = 0.9;

if alpha<0
    error('please set alpha between 0 and 1')
elseif alpha>1
    error('please set alpha between 0 and 1')
end

% position of source

p = 2;
q = 7;
r = 1;

if p>Lx
    error('please ensure position of source is within room')
elseif q>Ly
    error('please ensure position of source is within room')
elseif r>Lz
    error('please ensure position of source is within room')
end

% position of receiver

a = 2;
b = 4;
c = 5;

if a>Lx
    error('please ensure position of receiver is within room')
elseif b>Ly
    error('please ensure position of receiver is within room')
elseif c>Lz
    error('please ensure position of receiver is within room')
end

% IR calculation initializations

DELN = zeros((2*N+1),(2*N+1),(2*N+1));
G = DELN;
Nspan = (-N:N);
Aodd = Nspan*Lx-p-a;
Aeven = Nspan*Lx+p-a;
Bodd = Nspan*Ly-q-b;
Beven = Nspan*Ly+q-b;
Codd = Nspan*Lz-r-c;
Ceven = Nspan*Lz+r-c;
delnf = sqrt((((N+1)*Lx)^2)+(((N+1)*Ly)^2)+(((N+1)*Lz)^2))*Fs/343;
h = zeros(ceil(delnf),1);
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
            DELN(d+N+1,e+N+1,f+N+1) = round(t); % entire DELN array kept only for QC purposes
            G(d+N+1,e+N+1,f+N+1) = ((alpha)^w)/l; % entire G array kept only for QC purposes
            h(DELN(d+N+1,e+N+1,f+N+1)) = h(DELN(d+N+1,e+N+1,f+N+1)) + G(d+N+1,e+N+1,f+N+1);
        end
    end
end

% All reflections included check

if delnf<max(max(max(DELN)))
    disp ('DELN bad...somewhere')
elseif delnf>max(max(DELN))
    disp ('DELN all good.')
end

% plot DELN as a spatial QC. Animation moves through
% dimension varied by i

% for i=1:21
%     stem3(-N:N,-N:N,squeeze(DELN(i,:,:)))
%     drawnow
%     pause(0.1)
% end

% plot G as a spatial QC. Animation moves through
% dimension varied by i

% for i=1:21
%     surf(-N:N,-N:N,squeeze(G(i,:,:)))
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
alphastr = num2str(alpha);
Nstr = num2str(N);
params = strcat(Lxstr,'x',Lystr,'x',Lzstr,'_',alphastr,'_N',Nstr);

%output IR as 16-bit 44.1kHz wav

audiowrite(['ir_',params,'.wav'],h,Fs,'BitsPerSample',16);




