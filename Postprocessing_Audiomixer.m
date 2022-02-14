%Close existing graphs and clear namespace
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Adjust values below as desired

%Read in the song, other provided .mp3 files will work
mySong2 = audioread('pop.mp3');
info = audioinfo('pop.mp3');

%Experiment with different width of bands
%The width of each frequency band in Hz
fwidth = 200;

%Strengths array can be made different lengths. An error will throw if 
%   size of the totalsamples < array * fwidth
%Must also have at least 2 blocks

%Select one of the filters below or create your own

%All-pass
%strengths = [1,1];

%allstop
%strengths = [0,0];

%Lowpass (eliminate voices)
%strengths = [2,1.5,1,0.4,0,0];

%Equal response
%strengths = [0.04,0.25,1,1.5,2,2.5,3,5,1500,2000,2000];

%Random values
%strengths = [1.2,1.0,1.3,1.1,1.7,1.5,1.2,1,0.8,0.6,0.4];

%highpass (eliminate bass)
strengths = [0,0.2,0.4,1,1.5,2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
    % Program Run

%adjust information since reporting is different across both functions
info.TotalSamples = length(mySong2);
info.Duration = info.TotalSamples/info.SampleRate;

%Average the channels
%Assuming a small band and close together microphones this works
%When audio intentionally has noises on the different channels certain
%frequencies may get attenuated through averaging
mySong = zeros(info.TotalSamples, 1);
for i = 1:info.NumChannels
    mySong = mySong + mySong2(:,i);
end
mySong = mySong / info.NumChannels;

%Create time axis
t = 0:1/info.SampleRate:info.Duration;
t = t(1:end-1);

%Create f axis
f = info.SampleRate*(0:(info.TotalSamples/2))/info.TotalSamples;
f1 = info.SampleRate*(0:info.TotalSamples-1)/info.TotalSamples;

%Take the fourier transform of the signal
dfourier = fft(mySong);

%Scale specified frequencies to appropriate units
width = info.TotalSamples / info.SampleRate * fwidth;

%Compute filter
filter = zeros(info.TotalSamples, 1);
for i=1:size(strengths,2)-1
    filter(width*(i-1)+1:width*i+1) = strengths(i);
    filter(end-width*i:end-width*(i-1)) = strengths(i);
end
%Take last filter value to max w
    filter(width*((size(strengths,2))-1):end/2-1) = strengths(size(strengths,2));
    filter(end/2+1:end-width*(size(strengths,2)-1)) = strengths(size(strengths,2));

%Apply the filter to the frequency domain
newdfourier = dfourier.*filter;

%Reverse Fourier Transform to reconstitute the time-domain signal
ifft(newdfourier);
newsound = real(ifft(newdfourier));

%Play the filtered audio
sound(newsound, info.SampleRate)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting relavent graphs

%Plot multi-channel vs average-channel signals
figure
subplot(1,2,1)
plot(t, mySong2)
title('Multi-channel Volume')
xlabel('Time (s)')
ylabel('Audio Signal')
xlim([0,t(end)])
ylim([-1,1])

subplot(1,2,2)
plot(t, mySong)
title('Average-channel Volume')
xlabel('Time (s)')
ylabel('Audio Signal')
xlim([0,t(end)])
ylim([-1,1])
    
%Plot the filter
filterDisplay = filter(1:info.TotalSamples/2+1); %Remove negative frequencies
filterDisplay(2:end-1) = filterDisplay(2:end-1); %Adjust by one sample
figure
semilogx(f,filterDisplay) 
title('Filter')
xlabel('Frequency (Hz) [logarithmic axis]')
ylabel('Magnitude')
xlim([1,f(end)])

P2 = abs(dfourier/info.TotalSamples); %%Normalize
P1 = P2(1:info.TotalSamples/2+1); %%Remove negative frequencies
P1(2:end-1) = P1(2:end-1); %Adjust size to be accurate

%Plot the frequency domain vs the new frequency domain
figure
subplot(1,2,1)
semilogx(f,P1) 
title('Frequency Spectrum')
xlabel('Frequency (Hz) [logarithmic scale]')
ylabel('Relative Magnitude')
xlim([1,f(end)])

P2 = abs(newdfourier/info.TotalSamples); %%Normalize
P1 = P2(1:info.TotalSamples/2+1); %%Remove negative frequencies
P1(2:end-1) = P1(2:end-1); %Adjust size to be accurate
subplot(1,2,2)
semilogx(f, P1)
title('Filtered Frequency Spectrum')
xlabel('Frequency (Hz) [logarithmic scale]')
ylabel('Relative Magnitude')
xlim([1,f(end)])

%Plot original signal vs Filtered signal
figure
subplot(1,2,1)
plot(t, mySong)
title('Original Signal')
xlabel('Time (s)')
ylabel('Audio Signal')
xlim([0,t(end)])
ylim([-1,1])

subplot(1,2,2)
plot(t, newsound)
title('Transformed Signal')
xlabel('Time (s)')
ylabel('Audio Signal')
xlim([0,t(end)])
ylim([-1,1])