%% Some demonstrations of Signal Processing Fundamentals
%
% MDT
% 2016.01.12

%% Sine Waves
%
% The basic elements of signal analysis are sinusoidal waves. These are
% made with the sin() and cos() functions. Cosines are actually more
% fundamental mathematically, but as sines and cosines are the same (except
% for a shift) that is not a big deal. In the figure here, you can see that
% moving the cosine wave to the right (or sin to the left) will make the
% waves overlap.

timeLine = 0:0.001:2;    % From 0 to 2 seconds in steps 0f 0.001 (microseconds)
freq     = 1;            % in Hz (cps; cycles per second)

sw1 = sin(2*pi*freq*timeLine);
sw2 = cos(2*pi*freq*timeLine);

% Make a plot of a sin and a cos

plot(timeLine, sw1, 'b');
hold on;
plot(timeLine, sw2, 'r');
legend('Sine', 'Cosine', 'location', 'northeast');
title('Sine versus Cosine');
xlabel('Time (Seconds)');
ylabel('Amplitude (Arbitrary Units)');
grid on;
hold off;

%% Complex Waves
%
% The FFT "undoes" the wave from its complex form back into a breakdown of
% its component parts. There are a lot of tedious details, so don't lose
% the main point!

f1   = 0.50;        % Frequency of 1/2 Hz or 1/2 of a cycle per second (cps)
f2   = 2;           % 2 cps
tl   = 0:0.01:5;    % 5 seconds of time, measured in milliseconds

wave1 = sin(2*pi*f1*tl);
wave2 = sin(2*pi*f2*tl);
wave3 = wave1 + wave2;

% Whole waves:

plot(tl, wave1, 'b:');
hold on;
plot(tl, wave2, 'r--');
plot(tl, wave3, 'k');
legend('1 cycle per 2 sec', '2 cycles per sec', 'sum');
title('Sum of 2 waves');
xlabel('Time (Seconds)');
ylabel('Amplitude (Arbitrary Units)');
grid on;
hold off;

%% FFT
%
% The Fourier transform, usually implemented in terms of the "Fast" Fourier
% Transform (or FFT) takes a complex wave and reduces it to its component
% parts. Specifically, it determines the **frequencies** of the component
% sinusoid waves. (For instance, wave3 above would yield estimates of 1/2
% Hz and 2 Hz because those are the waves that go into making it.) 
%
% However, the FFT is for digitally sampled waves so there are some other
% details. For instance, the sampling rate (samples per second) and the
% length (duration) of the signal affect the results. This is why the
% calculations below are more complex.
%
% In the examples that follow, don't get lost in the details of making the
% plots -- the basic idea is to get the amplitude at each frequency. The
% problems include:
%
% * The peaks are not exactly at the frequencies we expect, but
%   instead their amplitude is spread out across nearby frequencies.
% * The LENGTH of the signal and the sampling rate define the amount
%   of smearing. (Note the sampling period especially.)

%% FFT Example 1 -- EEG Sampling Rate
%

Fs = 128;        % This is the sampling frequency (or sampling rate) for EEG
                 % 128 samples per second (or 128 Hz)
                 
T  = 1/Fs;       % Sampling period (length of a sample, 0.0078 sec or 7.8 ms)
L  = 1000;       % Length of signal (in samples)
t  = (0:L-1)*T;  % Just a trick for making a timeline of the right size (seconds)

% Make a complicated signal:

S = 1.0*sin(2*pi*2*t) + 1.0*sin(2*pi*10*t);  % 2 Hz and 10 Hz, equal amplitude

plot(t, S);
title('Complex Signal: 2 Hz + 10 Hz');
xlabel('Time (seconds)');
ylabel('S(t)');

% FFT

Y = fft(S);

% That was the easy part -- unfortunately we now have a vector of complex
% numbers! This leads to issues with making a figure. There are several
% standard ways to show the results, we use these:

% Two-sided Spectrum

P2 = abs(Y/L);                % Appropriately combine imaginary and real parts
P1 = P2(1:L/2+1);             % Take lower half of P2 (we want the One-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);  % Double all the values except for the first one

% One-sided Plot

f = Fs*(0:(L/2))/L;   % This is the frequency axis (look at it)
plot(f,P1)
title('Single-Sided Amplitude Spectrum of S(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% FFT Example 2 -- Higher Sampling Rate (from MATLAB Manual)
%
% This is a sample modified from the Matlab documentation on the fft
% function. It shows three simple waves of different frequencies and then
% shows the spectrum plot for each as a comparison.

Fs  = 1000;         % Sampling frequency (1000 samples per sec)
T   = 1/Fs;         % Sampling period (seconds per sample)
L   = 1000;         % Length of signal (samples)
t   = (0:L-1)*T;    % Time vector (seconds)

% Three waves at 50, 150, and 300 Hz:

x1 = cos(2*pi*50*t);    % First row wave
x2 = cos(2*pi*150*t);   % Second row wave
x3 = cos(2*pi*300*t);   % Third row wave
X  = [x1; x2; x3];      % Stick them into a container (for plotting)

% Plot the individual waves (for loop is just clever way to do it)

figure;
for i = 1:3
    subplot(3,1,i)
    plot(t(1:100),X(i,1:100))
    title(['Row ',num2str(i),' in the Time Domain'])
end

n = 2^nextpow2(L);  % WTF? FFT likes vectors of lengths equal to a power of 2.
                    %   Here we set the FFT length to 1024 = 2^10. 

% The weird bits below are because we stuck all of the waves into the 
% big-X container; basically that allows us to process three waves all at
% the same time rather than repeating all of this three times.

Y  = fft(X, n, 2);   % The "2" says each row in X is a different signal

P2            = abs(Y/n);           % Same as in previous section above
P1            = P2(:,1:n/2+1);
P1(:,2:end-1) = 2*P1(:,2:end-1);

% ...and make the one-sided spectrum plots:

figure;
for i=1:3
    subplot(3,1,i)
    plot(0:(Fs/n):(Fs/2-Fs/n),P1(i,1:n/2))
    title(['Row ',num2str(i), ' in the Frequency Domain'])
end

%% Noise
%
% Most signals have noise, such as motion artifacts, electrical lines, eye
% blinks, currents due to sweat, etc. in EEG. In this last example, we'll
% look at a relatively simple signal with some noise added. Here the noise
% will be plain old additive Gaussian noise (noise taken from a normal
% distribution). It is called additive because we will add it to the
% signal. :-)
%
% Note that the frequencies choen here are 25 and 100 Hz which will be easy
% to see in the Single-Sided Spectrum plot. Make sure this is completely
% clear when you are done.

close all; % This command just closes all the previous figures

%% Noise part I: The Signal (No Noise)
%

Fs = 1024;       % This is the INTERNAL sampling frequency for our EEG
T  = 1/Fs;       % Sampling period (length of a sample, 0.00098 s or 1 ms (rounded)
L  = 10000;      % Length of signal (in samples)
t  = (0:L-1)*T;  % A timeline of L samples (covers about 9.7 s)

% Make a complicated signal:

S = 5.0*sin(2*pi*25*t) + 10.0*sin(2*pi*100*t); % 25 & 100 Hz, Unequal amplitudes

% Hard to see, so we will look at only a small part of the wave:

figure;
plot(t(1:80), S(1:80));
title('Complex Signal (S): 25 and 100 Hz, Unequal Amplitudes');
xlabel('Time (seconds)');
ylabel('S(t)');

% FFT, Two-sided Spectrum, then One-sided Spectrum for Plotting

Y           = fft(S);         % FFT
P2          = abs(Y/L);       % Appropriately combine imaginary and real parts
P1          = P2(1:L/2+1);    % Take lower half of P2 (we want the One-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);  % Double all the values except for the first one
f           = Fs*(0:(L/2))/L; % This is the frequency axis (look at it)

% One-sided Plot

figure;
plot(f,P1)
title('Single-Sided Amplitude Spectrum of S(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% Noise part II: The Signal Plus the Noise
%
% Now we will add some noise to the signal. Literally.

noise_mu = 4;
noise_sd = 5;
noise    = noise_mu + noise_sd*randn(1,length(S));  % N(mu, sd)
SN       = S + noise;

% Redo everything we did above to get the plots of SN:

figure;
plot(t(1:80), SN(1:80));
title('Noisy Complex Signal (SN): 25 and 100 Hz, Unequal Amplitudes');
xlabel('Time (seconds)');
ylabel('SN(t)');

% FFT, Two-sided Spectrum, then One-sided Spectrum for Plotting

Y           = fft(SN);        % FFT
P2          = abs(Y/L);       % Appropriately combine imaginary and real parts
P1          = P2(1:L/2+1);    % Take lower half of P2 (we want the One-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);  % Double all the values except for the first one
f           = Fs*(0:(L/2))/L; % This is the frequency axis (look at it)

% One-sided Plot

figure;
plot(f,P1)
title('Single-Sided Amplitude Spectrum of SN(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% End
%
% This is just to close all the plots and clear Matlab out.

close all;
clear all;
