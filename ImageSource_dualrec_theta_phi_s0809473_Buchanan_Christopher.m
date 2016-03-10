clear all;

% this script approximates presence of two ears (or receivers), outputting an IR for
% each - imaginatively labelled with an L and R. Twiddling the r, theta, and phi knobs resizes and rotates
% head in xy, xz and yz planes.

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

% position of receivers' mid-point

a = 2;
b = 4;
c = 5;

% radius (r) and orientation (theta) of head in xyz space.

r = 0.5;
theta = pi/2;
phi = pi/2;

if abs(theta)>2*pi
    error('theta within 2*pi please')
elseif abs(phi)>2*pi
    error('phi within 2*pi please')
end

% parametric solutions of a sphere in 3D Cartesian

al = r*cos(theta)*sin(phi)+a;
ar = r*cos(theta+pi)*sin(phi)+a;
bl = r*sin(theta)*sin(phi)+b;
br = r*sin(theta+pi)*sin(phi)+b;
cl = r*cos(phi)+c;
cr = r*cos(phi+pi)+c;

if al>Lx
    error('please ensure position of receivers are within room')
elseif ar>Lx
    error('please ensure position of receivers are within room')
elseif bl>Ly
    error('please ensure position of receivers are within room')
elseif br>Ly
    error('please ensure position of receivers are within room')
elseif cl>Lz
    error('please ensure position of receivers are within room')
elseif cr>Lz
    error('please ensure position of receivers are within room')
end


% IR calculation initializations

DELNl = zeros((2*N+1),(2*N+1),(2*N+1));
DELNr = DELNl;
Gl = DELNl;
Gr = Gl;
Nspan = (-N:N);
Aoddl = Nspan*Lx-p-al;
Aoddr = Nspan*Lx-p-ar;
Aevenl = Nspan*Lx+p-al;
Aevenr = Nspan*Lx+p-ar;
Boddl = Nspan*Ly-q-bl;
Boddr = Nspan*Ly-q-br;
Bevenl = Nspan*Ly+q-bl;
Bevenr = Nspan*Ly+q-br;
Coddl = Nspan*Lz-r-cl;
Coddr = Nspan*Lz-r-cr;
Cevenl = Nspan*Lz+r-cl;
Cevenr = Nspan*Lz+r-cr;
delnf = sqrt((((N+1)*Lx)^2)+(((N+1)*Ly)^2)+(((N+1)*Lz)^2))*Fs/343;
hl = zeros(ceil(delnf),1);
hr = hl;
ti = Fs/343;

% DELN (delay times in samples) and G (amplitudes) computation

for d=-N:N;
    for e=-N:N;
        for f=-N:N;
            if mod(d,2)==1
                Al = Aoddl(d+N+2);
                Ar = Aoddr(d+N+2);
            else
                Al = Aevenl(d+N+1);
                Ar = Aevenr(d+N+1);
            end
            if mod(e,2)==1
                Bl = Boddl(e+N+2);
                Br = Boddr(e+N+2);
            else
                Bl = Bevenl(e+N+1);
                Br = Bevenr(e+N+1);
            end
            if mod(f,2)==1
                Cl = Coddl(f+N+2);
                Cr = Coddr(f+N+2);
            else
                Cl = Cevenl(f+N+1);
                Cr = Cevenr(f+N+1);
            end
            ll = sqrt((Al^2)+(Bl^2)+(Cl^2));
            lr = sqrt((Ar^2)+(Br^2)+(Cr^2));
            tl = ti*ll;
            tr = ti*lr;
            w = abs(d)+abs(e)+abs(f);
            DELNl(d+N+1,e+N+1,f+N+1) = round(tl); % entire DELN array kept only for QC purposes
            DELNr(d+N+1,e+N+1,f+N+1) = round(tr);
            Gl(d+N+1,e+N+1,f+N+1) = ((alpha)^w)/ll; % entire G array kept only for QC purposes
            Gr(d+N+1,e+N+1,f+N+1) = ((alpha)^w)/lr;
            hl(DELNl(d+N+1,e+N+1,f+N+1)) = hl(DELNl(d+N+1,e+N+1,f+N+1)) + Gl(d+N+1,e+N+1,f+N+1);
            hr(DELNr(d+N+1,e+N+1,f+N+1)) = hr(DELNr(d+N+1,e+N+1,f+N+1)) + Gr(d+N+1,e+N+1,f+N+1);
        end
    end
end

% QC of head size and orientation relative to room

xroom = [0,0,0,0,Lx,Lx,Lx,Lx];
yroom = [0,Ly,0,Ly,Ly,0,Ly,0];
zroom = [0,0,Lz,Lz,Lz,Lz,0,0];
plot3(xroom,yroom,zroom,'r');
grid on
hold on
xrec = [al,ar];
yrec = [bl,br];
zrec = [cl,cr];
scatter3(xrec,yrec,zrec,'filled');

% All reflections included QC

if delnf<max(max(max(DELNl)))
    disp ('DELNl bad...somewhere')
elseif delnf>max(max(DELNl))
    disp ('DELNl all good.')
end
if delnf<max(max(max(DELNr)))
    disp ('DELNr bad...somewhere')
elseif delnf>max(max(DELNr))
    disp ('DELNr all good.')
end

% plot DELNl as a spatial QC. Animation moves through
% dimension varied by i

% for i=1:21
%     stem3(-N:N,-N:N,squeeze(DELNl(i,:,:)))
%     drawnow
%     pause(0.1)
% end

% plot DELNr as a spatial QC. Animation moves through
% dimension varied by i

% for i=1:21
%     stem3(-N:N,-N:N,squeeze(DELNr(i,:,:)))
%     drawnow
%     pause(0.1)
% end

% plot Gl as a spatial QC. Animation moves through
% dimension varied by i

% for i=1:21
%     surf(-N:N,-N:N,squeeze(Gl(i,:,:)))
%     drawnow
%     pause(0.1)
% end

% normalize IRs

hl = hl/(norm(hl,Inf));
ll = length(hl);
hr = hr/(norm(hr,Inf));
lr = length(hr);

% plot impulse responses

% T = (0:Ts:(ll-1)*Ts)';
% plot(T,hl,T,hr)
% xlabel('time (s)')
% ylabel('amplitude')
% title('IR')

%generate impulse response filenames

Lxstr = num2str(Lx);
Lystr = num2str(Ly);
Lzstr = num2str(Lz);
alphastr = num2str(alpha);
Nstr = num2str(N);
params = strcat(Lxstr,'x',Lystr,'x',Lzstr,'_',alphastr,'_N',Nstr);

%output IR as 16-bit 44.1kHz wav

audiowrite(['ir_',params,'_l.wav'],hl,Fs,'BitsPerSample',16);
audiowrite(['ir_',params,'_r.wav'],hr,Fs,'BitsPerSample',16);




