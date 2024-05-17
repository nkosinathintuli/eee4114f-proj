% Import a WAV file
[originalAudio, fs] = audioread('harvard.wav');

% Play the original audio
disp('Playing original audio...');
%sound(originalAudio, fs);
%pause(length(originalAudio)/fs + 1);

% Add white noise to the audio
noiseLevel = 0.05; % Adjust this value for different noise levels
noisyAudio = originalAudio + noiseLevel * randn(size(originalAudio));

% Play the noisy audio
disp('Playing noisy audio...');
%sound(noisyAudio, fs);
%pause(length(noisyAudio)/fs + 1);

% Design an IIR Butterworth filter to remove the white noise
filterOrder = 6; % Order of the Butterworth filter
cutoffFreq = 2000; % inspection of the audio 
cutoffFreqNorm = 2 * cutoffFreq/fs; % Normalized cutoff frequency (0 to 1, where 1 is Nyquist frequency)
[b, a] = butter(filterOrder, cutoffFreqNorm, 'low');


% Apply the Butterworth filter to the noisy audio
filteredAudio = filter(b, a, noisyAudio);

% Play the filtered audio
disp('Playing filtered (denoised) audio...');
%sound(filteredAudio, fs);
%pause(length(filteredAudio)/fs + 1);

freqz(b,1);

% Plot the original, noisy, and filtered audio signals for comparison
t = (0:length(originalAudio)-1) / fs;
figure;
subplot(3,1,1);
plot(t, originalAudio);
title('Original Audio');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3,1,2);
plot(t, noisyAudio);
title('Noisy Audio');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3,1,3);
plot(t, filteredAudio);
title('Filtered Audio');
xlabel('Time (s)');
ylabel('Amplitude');


% Compute and display SNR, PSNR, MSE values
mse_noisy = mean((originalAudio - noisyAudio).^2);
mse_filtered = mean((originalAudio - filteredAudio).^2);
mse_improvement = ((mse_noisy - mse_filtered) / abs(mse_noisy)) * 100;
mse_noisy=mse_noisy(1);
mse_filtered=mse_filtered(1);


peakSignal = max(abs(originalAudio));


snr_noisy = sum(originalAudio.^2) / sum((originalAudio - noisyAudio).^2)
snr_noisy_dB = 10 * log10(snr_noisy);

snr_filtered = sum(originalAudio.^2) / sum((originalAudio - filteredAudio).^2)
snr_filtered_dB = 10 * log10(snr_filtered);

snr_improvement = ((snr_filtered - snr_noisy) / abs(snr_noisy)) * 100;

psnr_noisy = 20 * log10(peakSignal / sqrt(mse_noisy));
psnr_filtered = 20 * log10(peakSignal / sqrt(mse_filtered));
psnr_improvement = ((psnr_filtered - psnr_noisy) / abs(psnr_noisy)) * 100;
psnr_noisy=psnr_noisy(1);
psnr_filtered=psnr_filtered(1);

fprintf("\n");
fprintf('MSE (Noisy): %f\n', mse_noisy);
fprintf('MSE (Filtered): %f\n', mse_filtered);
fprintf('Percentage Improvement in MSE: %.2f%%\n\n', mse_improvement);

fprintf('SNR (Noisy): %f dB\n', snr_noisy_dB);
fprintf('SNR (Filtered): %f dB\n', snr_filtered_dB);
fprintf('Percentage Improvement in SNR: %.2f%%\n\n', snr_improvement);

fprintf('PSNR (Noisy): %f dB\n', psnr_noisy);
fprintf('PSNR (Filtered): %f dB\n', psnr_filtered);
fprintf('Percentage Improvement in PSNR: %.2f%%\n', psnr_improvement);
