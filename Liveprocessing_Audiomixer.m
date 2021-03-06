%Clear old stuff
clear all;
close all;

%Multiplier for each filter block
bcoeffs = [1,1,1,1,1,1];

% NOTE
% Due to the way the code is implemented, for the effects of changing the
% values for the different frequencies to take place in the program, the
% pause button must be pressed, and then the resume button. You don't have
% to wait at all, so you could just double click the button. It may take a
% second or 2 for the new effects to buffer to the song.
% NOTE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GUI

%Create the window
S.fh = figure('units','pixels',...
              'position',[100 100 200 320],...
              'menubar','none',...
              'name','GUI_15',...
              'numbertitle','off',...
              'resize','off');
      
%Create all of the text boxes
%Lowest box
S.ed1 = uicontrol('style','edit',...
                 'unit','pix',...
                 'position',[90 70 80 30],...
                 'string','1.0');
S.tx1 = uicontrol('style','text',...
                 'unit','pix',...
                 'position',[10 70 80 30],...
                 'string','<500 Hz',...
                 'fontsize',12);
             
S.ed2 = uicontrol('style','edit',...
                 'unit','pix',...
                 'position',[90 110 80 30],...
                 'string','1.0');
S.tx2 = uicontrol('style','text',...
                 'unit','pix',...
                 'position',[10 110 80 30],...
                 'string','0.5-1 Hz',...
                 'fontsize',12);
             
S.ed3 = uicontrol('style','edit',...
                 'unit','pix',...
                 'position',[90 150 80 30],...
                 'string','1.0');
S.tx3 = uicontrol('style','text',...
                 'unit','pix',...
                 'position',[10 150 80 30],...
                 'string','1-1.5 kHz',...
                 'fontsize',12);
             
S.ed4 = uicontrol('style','edit',...
                 'unit','pix',...
                 'position',[90 190 80 30],...
                 'string','1.0');
S.tx4 = uicontrol('style','text',...
                 'unit','pix',...
                 'position',[10 190 80 30],...
                 'string','1.5-2 kHz',...
                 'fontsize',12);
             
S.ed5 = uicontrol('style','edit',...
             'unit','pix',...
             'position',[90 230 80 30],...
             'string','1.0');
S.tx5 = uicontrol('style','text',...
                 'unit','pix',...
                 'position',[10 230 80 30],...
                 'string','2-2.5 kHz',...
                 'fontsize',12);
             
S.ed6 = uicontrol('style','edit',...
             'unit','pix',...
             'position',[90 270 80 30],...
             'string','1.0');
S.tx6 = uicontrol('style','text',...
                 'unit','pix',...
                 'position',[10 270 80 30],...
                 'string','2.5< kHz',...
                 'fontsize',12);

%Create the pushbutton
S.pb = uicontrol('style','push',...
                 'unit','pix',...
                 'position',[10 30 180 30],...
                 'Interruptible','off',...
                 'string','Button Does Nothing'); 
%Set up callback function for the button
set([S.pb],'callback',{@pb_call,S}) % Set callbacks.       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Wait a second to ensure the GUI gets set up before it is locked out
pause(1);

%Load in preconstructed bandpass filter blocks
b1 = load('block1.mat');
b2 = load('block2.mat');
b3 = load('block3.mat');
b4 = load('block4.mat');
b5 = load('block5.mat');
b6 = load('block6.mat');

%Add coefficients of the individual difference equations to get a combined
%filter
b = bcoeffs(1)*b1.b1.Coefficients + bcoeffs(2)*b2.b2.Coefficients + bcoeffs(3)*b3.b3.Coefficients + bcoeffs(4)*b4.b4.Coefficients + bcoeffs(5)*b5.b5.Coefficients + bcoeffs(6)*b6.b6.Coefficients;

%insert link to website and this line might work, needs to be a link
%directly to a audio file
%stream = dsp.AudioFileReader('http://audio.wgbh.org:8004/')

%SamplesPerFrame means the system sends a bunch of samples at once. For
%smaller frames the audio quality is greatly decreased for some reason
stream = dsp.AudioFileReader('chill.mp3', 'SamplesPerFrame', 32768);

%instantiate writer device to the audio drivers
writer = audioDeviceWriter('SampleRate', stream.SampleRate);

%Continue outputting until the stream sends 1 instead of 0 to signify
%termination
while ~isDone(stream)
    %Check to see if the user has adjusted the frequency bands
    bcoeffs(1) = str2double(get(S.ed1, 'string'));
    bcoeffs(2) = str2double(get(S.ed2, 'string'));
    bcoeffs(3) = str2double(get(S.ed3, 'string'));
    bcoeffs(4) = str2double(get(S.ed4, 'string'));
    bcoeffs(5) = str2double(get(S.ed5, 'string'));
    bcoeffs(6) = str2double(get(S.ed6, 'string'));
    %Update the difference equation accordingly
    b = bcoeffs(1)*b1.b1.Coefficients + bcoeffs(2)*b2.b2.Coefficients + bcoeffs(3)*b3.b3.Coefficients + bcoeffs(4)*b4.b4.Coefficients + bcoeffs(5)*b5.b5.Coefficients + bcoeffs(6)*b6.b6.Coefficients;
    
    %Get the audio input as a 32768x2 array    
    audio = stream();
    
    %Filter the input with the difference equation, y[n] = ax[n] + bx[n-1]
    % + ... + zx[n-442]
    y3 = filter(b,1,audio);
    
    %Write the new output back as a large frame to the drivers. If the
    %drivers buffer size is too small, issues may arise, large buffers
    %cause catch up latency between the pause button and the music stopping
    writer(y3);
end

%Release locks from computer
release(stream);
release(writer);


%Callback function executes when the button is pressed. Has lower priority
%than the audio stream so only works when the program is paused
function [] = pb_call(varargin)
S = varargin{3};
%disp(get(S.ed1,'string')) 
end

