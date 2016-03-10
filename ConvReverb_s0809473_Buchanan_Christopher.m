clear all

% Have user browse for a file, starting at current location.

currentFolder = pwd;

% Get the full path of the audio file that the user wants to reverberate.

audiofilepath = fullfile(currentFolder, '*.*');
disp('Pick an audio sample please')
[audiofilename, folder] = uigetfile(audiofilepath, 'Pick a sample to reverberate');
if audiofilename == 0  % user clicked cancel (maybe because of stage fright).
    error('Fine then. Dont use my script, see if I care.');
end
filepathin = fullfile(folder, audiofilename);

[x,Fs] = audioread(filepathin);

if size(x,2)>1
    x = (x(:,1)+x(:,2))/2;
    disp('input is stereo - converted to mono')
end

Nyq = Fs/2;

% read in IR

IRfilepath = fullfile(currentFolder, '*.*');
disp('Pick an IR please')
[IRfilename, folder] = uigetfile(IRfilepath, 'Pick an IR');
if IRfilename == 0
   error('Fine then. Dont use my script, see if I care.');
end
filepathin = fullfile(folder, IRfilename);
[h,Fsh] = audioread(filepathin);

[~,audio,aext] = fileparts(audiofilename);
[~,IR,irext] = fileparts(IRfilename);

% compatibility checks

if aext~=irext
    error('file formats do not match. Not sure how to deal with this yet.')
end

if size(h,2)>1
    error('This is a stereo IR. Not sure how to deal with this yet.')
end

if Fs ~= Fsh
    error('sample rates do not match. Not sure how to deal with this yet.')
end

% time domain convolution (warning: SLOW)

% l= length(x);       
% k = length(h); 
% xpad = [x;zeros(k,1)];  
% hpad = [h;zeros(l,1)];  % zero padding of both x[n] and h[n] if size of x[n] ~= size of h[n]
% y = zeros(length(xpad),1);
% for j = 1:l+k-1 
%      for i = 1:l 
%          if(j-i+1>0) 
%              y(j) = y(j) + xpad(i)*hpad(j-i+1);
%          else 
%          end 
%      end
% end 

% frequency domain convolution

l= length(x);       
k = length(h); 
xpad = [x;zeros(k,1)];  
hpad = [h;zeros(l,1)];  % zero padding of both x[n] and h[n] if size of x[n] ~= size of h[n]
XF = fft(xpad);
HF = fft(hpad);
YF = XF.*HF;
y = ifft(YF);
y = y/norm(y,Inf);

% convolution using conv

% y = conv(x,h);
% y = y/norm(y,Inf);

% output

audiowrite([audio,'_',IR,'.wav'],y,Fs);
sound(y,Fs);




