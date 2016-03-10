clear all;

% This script adds mono reverb to an audio file of the user's choosing,
% outputting it with a '_reverbed_ex' suffix. Any number of FFC filters can
% be implemented here, thus in general "b" and "M" for the feedforcombex
% function are row vectors, with a straightforward 1-1 correspondence
% between entries in b and entries in M.

% Have user browse for a file, starting at current location.
currentFolder = pwd;

% Get the full path of the audio file that the user wants to Schroeder.
defaultfilename = fullfile(currentFolder, '*.*');
[filename, folder] = uigetfile(defaultfilename, 'Pick something');
if filename == 0  % user clicked cancel (maybe because of stage fright).
    return;
end
filepathin = fullfile(folder, filename);

[x,Fs] = audioread(filepathin);

if size(x,2)>1
    x = (x(:,1)+x(:,2))/2;
    disp('input is stereo - converted to mono')
end

Nyq = Fs/2;

%AP1

y1 = allpassM_s0809473_Buchanan_Christopher(x,0.7,1051);

%AP2

y2 = allpassM_s0809473_Buchanan_Christopher(y1,0.7,337);

%AP3

y3 = allpassM_s0809473_Buchanan_Christopher(y2,0.7,113);

%FFCF1 - FFCF4

y4 = feedforcombex_s0809473_Buchanan_Christopher(y3,[0.742,0.733,0.715,0.697],[4799,4999,5399,5801]);

%normalize and output

yfinal3 = y4/norm(y4,Inf);

[~,name,ext] = fileparts(filename);

audiowrite([name,'_reverbed_ex',ext],yfinal3,Fs);

