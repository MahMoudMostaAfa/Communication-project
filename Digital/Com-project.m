%% QUESTION 1
% Define the time vector
% Define the sampling frequency and resolution
fs = 100;   % Sampling frequency (Hz)
df = 0.01;  % Resolution (Hz)

% Define the time vector
t = -2:1/fs:2;
N = length(t); % Length of the signal

% Define the signal x(t)
x = zeros(size(t));
x(t >= -2 & t <= -1) = t(t >= -2 & t <= -1) + 2;
x(t >= -1 & t <= 1) = 1;
x(t >= 1 & t <= 2) = -t(t >= 1 & t <= 2) + 2;

% Calculate the Fourier transform analytically
f = -fs/2:df:fs/2;  % Frequency vector
X_analytical = zeros(size(f));
X_analytical(abs(f) <= 1) = 1./(2 * pi * 1i * f(abs(f) <= 1));
X_analytical(f == 0) = 1;

% Calculate the Fourier transform numerically using FFT
X_numerical = fftshift(fft(x)) / N;

% Plot the signals and their Fourier transforms
ظظfigure;

subplot(2, 1, 1);
plot(t, x, 'b', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Signal x(t)');
grid on;

subplot(2, 1, 2);
plot(f, abs(X_analytical), 'r', 'LineWidth', 2);
hold on;
stem(f, abs(X_numerical), 'b', 'LineWidth', 1);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Fourier Transform');
legend('Analytical', 'Numerical');
grid on;
% Calculate the power spectrum of the signal
Pxx = abs(X_numerical).^2;

% Find the maximum power
P_max = max(Pxx);

% Find the frequency band where power drops to 5% of the maximum
BW_index = find(Pxx >= 0.05 * P_max, 1, 'last');
BW = 2 * abs(f(BW_index)); % BW in Hz
disp(['Estimated Bandwidth (BW): ', num2str(BW), ' Hz']);
% Define the cutoff frequencies for LPF
BW_1 = 1;   % Bandwidth of the LPF in Hz
BW_2 = 0.3; % Reduced bandwidth of the LPF in Hz

% Design LPF transfer functions
[b_1, a_1] = butter(6, BW_1/(fs/2), 'low'); % Butterworth LPF with BW_1
[b_2, a_2] = butter(6, BW_2/(fs/2), 'low'); % Butterworth LPF with BW_2

% Filter the signal using LPF
y_1 = filter(b_1, a_1, x);
y_2 = filter(b_2, a_2, x);

% Plot the filtered signals
figure;
subplot(2, 1, 1);
plot(t, x, 'b', 'LineWidth', 2);
hold on;
plot(t, y_1, 'r', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Input Signal and Filter Output (BW = 1 Hz)');
legend('Input Signal', 'Filtered Signal');
grid on;

subplot(2, 1, 2);
plot(t, x, 'b', 'LineWidth', 2);
hold on;
plot(t, y_2, 'r', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Input Signal and Filter Output (BW = 0.3 Hz)');
legend('Input Signal', 'Filtered Signal');
grid on;






