% MATLAB script to analyze an audio signal from a WAV file

% Clear workspace and close all figures
clear;
close all;
clc;

% Read the audio file
[audioSignal, sampleRate] = audioread('harvard.wav');
noiseLevel = 0.05;
noise = noiseLevel * randn(size(audioSignal));

% Get the length of the audio signal
n = length(audioSignal);


% Time vector for the audio signal
t = (0:n-1) / sampleRate;

% Plot the audio signal in time domain
figure;
plot(t, audioSignal);
xlabel('Time (s)');
ylabel('Amplitude');
title('Audio Signal in Time Domain');
grid on;

% Compute the Fourier Transform of the audio signal
Y = fft(audioSignal);
Noise = fft(noise);

% Compute the two-sided spectrum
P2 = abs(Y/n);
%P2 = abs(Noise/n);

% Compute the single-sided spectrum based on the even-valued signal length
P1 = P2(1:n/2+1);

P1(2:end-1) = 2*P1(2:end-1);

% Frequency vector
f = sampleRate*(0:(n/2))/n;

% Plot the single-sided amplitude spectrum
figure;
plot(f, P1);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Single-Sided Amplitude Spectrum of Audio Signal');
grid on;

% Optional: Zoom in on a specific frequency range (e.g., 0 to 5000 Hz)
%xlim([0 5000]);

saveas(gcf, 'FrequencySpectrum.png');