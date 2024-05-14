pkg load signal
clc;
clear all;
close all;
% Generate random binary data
m = randi([0 1], 1, 64); % 64-bit random binary data
n = length(m);
%Implementation of UNIPOLAR NRZ
% Parameters
max_amplitude = 1; % Maximum amplitude of the signal
min_amplitude = 0; % Minimum amplitude of the signal
T = 1000; % Total time of the signal
fs = 1000; % Sampling rate
ts = 1/fs; % Sampling time
t = linspace(0, n, n*fs); % Time vector

% Generate pulse train
y = [];
for i = 1:n
    if m(i) == 1
        y = [y linspace(max_amplitude, max_amplitude, fs)]; % Fill y with max values for this interval
    else
        y = [y linspace(min_amplitude, min_amplitude, fs)]; % Fill y with min values for this interval
    end
end

% Plot temporal characteristics of the transmitted signal
figure;
plot(t, y);
axis([0 n -2 2]); %make y axis plot from -2 to 2
grid on;
box off;
ylabel('Signal(V)');
xlabel('Time(s)');
title('Unipolar NRZ');

% ASK Modulation
fc =3; % Carrier frequency (Hz) higher than Rb , Rb=1
carrier = cos(2*pi*fc*t); % Carrier signal

% ASK Modulated Signal
ask_signal = y .* carrier;

% Plot ASK Modulated Signal
figure;
plot(t, ask_signal);
axis([0 n -3 3]);
grid on;
box off;
xlabel('Time (s)');
ylabel('Amplitude');
title('ASK Modulated Signal');

% ASK Receiver Implementation

% Define carrier frequency and threshold
threshold = 0.5; % Threshold for detection

% Generate local carrier signals with different phase shifts
carrier_30 = cos(2*pi*fc*t + deg2rad(30)); % 30 degree phase shift
carrier_60 = cos(2*pi*fc*t + deg2rad(60)); % 60 degree phase shift
carrier_90 = cos(2*pi*fc*t + deg2rad(90)); % 90 degree phase shift

% Coherent detection: multiply received signal with carriers having different phase shifts
received_signal_30 = ask_signal .* carrier_30; % Received ASK modulated signal with 30 degree phase shift
received_signal_60 = ask_signal .* carrier_60; % Received ASK modulated signal with 60 degree phase shift
received_signal_90 = ask_signal .* carrier_90; % Received ASK modulated signal with 90 degree phase shift

% Threshold detection
demodulated_data_30 = received_signal_30 > threshold;
demodulated_data_60 = received_signal_60 > threshold;
demodulated_data_90 = received_signal_90 > threshold;

% Low-pass filtering in the time domain
cutoff_frequency = 2; % Cutoff frequency for LPF (Hz)
[b, a] = butter(7, cutoff_frequency/(fs/2), 'low'); % Butterworth LPF coefficients

% Apply LPF to the demodulated signals
filtered_data_30 = filter(b, a, demodulated_data_30);
filtered_data_60 = filter(b, a, demodulated_data_60);
filtered_data_90 = filter(b, a, demodulated_data_90);

% Plotting
figure;

% Plot received ASK modulated signals
subplot(4,2,1);
plot(t, received_signal_30);
xlabel('Time (s)');
ylabel('Amplitude');
title('Received ASK Modulated Signal with 30° Phase Shift');

subplot(4,2,3);
plot(t, received_signal_60);
xlabel('Time (s)');
ylabel('Amplitude');
title('Received ASK Modulated Signal with 60° Phase Shift');

subplot(4,2,5);
plot(t, received_signal_90);
xlabel('Time (s)');
ylabel('Amplitude');
title('Received ASK Modulated Signal with 90° Phase Shift');

% Plot demodulated data after threshold detection
subplot(4,2,2);
plot(t, demodulated_data_30, 'r');
xlabel('Time (s)');
ylabel('Demodulated Data');
title('Demodulated Data (After Threshold Detection) with 30° Phase Shift');

subplot(4,2,4);
plot(t, demodulated_data_60, 'r');
xlabel('Time (s)');
ylabel('Demodulated Data');
title('Demodulated Data (After Threshold Detection) with 60° Phase Shift');

subplot(4,2,6);
plot(t, demodulated_data_90, 'r');
xlabel('Time (s)');
ylabel('Demodulated Data');
title('Demodulated Data (After Threshold Detection) with 90° Phase Shift');

% Plot filtered demodulated data
% Plot all filtered demodulated data in one figure
figure;

% Plot filtered demodulated data for 30° phase shift
subplot(3,1,1);
plot(t, filtered_data_30, 'g');
xlabel('Time (s)');
ylabel('Amplitude');
title('Filtered Demodulated Data with 30° Phase Shift');

% Plot filtered demodulated data for 60° phase shift
subplot(3,1,2);
plot(t, filtered_data_60, 'g');
xlabel('Time (s)');
ylabel('Amplitude');
title('Filtered Demodulated Data with 60° Phase Shift');

% Plot filtered demodulated data for 90° phase shift
subplot(3,1,3);
plot(t, filtered_data_90, 'g');
xlabel('Time (s)');
ylabel('Amplitude');
title('Filtered Demodulated Data with 90° Phase Shift');

% Adjust subplot spacing

% Fourier transform of each spectrum
figure;


% Compute spectrum
spectrum_orig = fft(ask_signal);
spectrum_30 = fft(filtered_data_30);
spectrum_60 = fft(filtered_data_60);
spectrum_90 = fft(filtered_data_90);
f = linspace(-fs/2, fs/2, length(filtered_data_30));

% Plot spectrum
subplot(4,1,1);
plot(f, fftshift(abs(spectrum_orig))*ts);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Fourier Transform of Original ASK Modulated Signal');

% Plot spectrum for filtered demodulated data with 30° phase shift
subplot(4,1,2);
plot(f, fftshift(abs(spectrum_30))*ts);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Fourier Transform of Filtered Demodulated Data with 30° Phase Shift');

% Plot spectrum for filtered demodulated data with 60° phase shift
subplot(4,1,3);
plot(f, fftshift(abs(spectrum_60))*ts);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Fourier Transform of Filtered Demodulated Data with 60° Phase Shift');

% Plot spectrum for filtered demodulated data with 90° phase shift
subplot(4,1,4);
plot(f, fftshift(abs(spectrum_90))*ts);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Fourier Transform of Filtered Demodulated Data with 90° Phase Shift');

