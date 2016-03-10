clear all

% this script is to be used with '_l' and '_r' outputs of
% ImageSource_dualrec scripts, to approximate experience of dual-receiver
% listener.

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

% read in IRs

IRlfilepath = fullfile(currentFolder, '*.*');
disp('Pick a LEFT IR please')
[IRlfilename, folder] = uigetfile(IRlfilepath, 'Pick a LEFT IR');
if IRlfilename == 0
   error('Fine then. Dont use my script, see if I care.');
end
filepathin = fullfile(folder, IRlfilename);
[hl,Fshl] = audioread(filepathin);

[~,audio,aext] = fileparts(audiofilename);
[~,IRl,irlext] = fileparts(IRlfilename);

IRrfilepath = fullfile(currentFolder, '*.*');
disp('Pick a RIGHT IR please')
[IRrfilename, folder] = uigetfile(IRrfilepath, 'Pick a RIGHT IR');
if IRrfilename == 0
   error('Fine then. Dont use my script, see if I care.');
end
filepathin = fullfile(folder, IRrfilename);

[hr,Fshr] = audioread(filepathin);

[~,IRr,irrext] = fileparts(IRrfilename);

% compatibility checks

if aext~=irlext
    error('file formats do not match. Not sure how to deal with this yet.')
end

if size(hl,2)>1
    error('This is a stereo IR. Not sure how to deal with this yet.')
end

if Fs ~= Fshl
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
kl = length(hl); 
xpad = [x;zeros(kl,1)];  
hlpad = [hl;zeros(l,1)];  
XF = fft(xpad);
HLF = fft(hlpad);
YLF = XF.*HLF;
yl = ifft(YLF);
yl = yl/norm(yl,Inf);
      
kr = length(hr); 
xpad = [x;zeros(kr,1)];  
hrpad = [hr;zeros(l,1)];  
XF = fft(xpad);
HRF = fft(hrpad);
YRF = XF.*HRF;
yr = ifft(YRF);
yr = yr/norm(yr,Inf);

% output

y = [yl,yr];
audiowrite([audio,'_',IRl,'_r','.wav'],y,Fs);
sound(y,Fs);




