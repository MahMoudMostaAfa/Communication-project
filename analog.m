clc;
close all
clear
pkg load signal;
pkg load symbolic;

%Initialiazing time and frequency Vectors%
fs = 100; % samples / sec
T = 100; % df = fs/N;
df = 1/T; %frequency step
N = T*fs;%N/T = 100
ts = 1/fs; %sampling time
t = linspace(-50,50,N); %time vector
if(rem(N,2)==0) %% Even
  f = - (0.5*fs) : df : (0.5*fs-df) ; %% Frequency vector if x/f is even
else %% Odd
  f = - (0.5*fs-0.5*df) : df : (0.5*fs-0.5*df) ; %% Frequency vector if x/f is odd
end

% Step_1 Plot the function x(t)%
x_t = zeros(1,N);
x_t(t>=-1 & t<=1) = 1;
x_t(t<=-1 & t>=-2)= t(t<=-1 & t>=-2) + 2;
x_t(t>=1 & t<=2) = 2 - t(t>=1 & t<=2);
% Plot the signal
figure;
plot(t, x_t),axis([-3 3 -2 2]);
title('Signal x(t)');
xlabel('t');
ylabel('x(t)');
grid on

% Step_2 Derive an analytical expression for its Fourier transform%
%FT(rect(pulsewidth = 1)).FT(rect(pulsewidth = 3))
X_f_analytic = 3*sinc(3*f).*sinc(f);

% Step_3 Plot analytical expression and Spectrum of sampled signal%
X_f = (fftshift(fft(x_t)))*ts ; % dt non-periodic signal
figure;
plot(f,abs(X_f_analytic),'-;analytic expression;','LineWidth',2, f,abs(X_f),'-;Spectrum of sampled signal;','LineWidth',2), axis([-3 3]);
title('X(f)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on
box off

% Step_4 Estimate the BW @95% power %
Power_max = max(abs(X_f).^2); %% Or  sum(abs(m).^2)*ts;
Index = find(f == 0);%% Index of zero frequency
for c_index = length(f):-1:length(f)/2
Power = abs(X_f(c_index).^2);
  if(Power >= 0.05*0.5*Power_max)
    BW_x = f(c_index);
    break
  end
end

% Step_5 Ideal LPF with f(cut-off) = 1Hz %
H1 = abs(f) <= 1;
Y1_f = H1 .* X_f;
Y1_t = ifft(ifftshift(Y1_f)/ts);
figure;
plot(t,Y1_t,'-;Filter Output y1(t);','LineWidth',2, t,x_t,'-;Filter Input;','LineWidth',2),axis([-5 5]);
title('x(t)');
xlabel('Time(sec)');
ylabel('Magnitude');
grid on
box off

% Step_6 Ideal LPF with f(cut-off) = 0.3Hz %
H2 = abs(f) <= 0.3;
Y2_f = H2.*X_f;
Y2_t = ifft(ifftshift(Y2_f)/ts);
figure;
plot(t,Y2_t,'-;Y2(t);','LineWidth',2, t,x_t,'-;Filter Input;','LineWidth',2),axis([-10 10]);
title('x(t)');
xlabel('Frequency(Hz)');
ylabel('Magnitude');
grid on
box off

% Step_7.1 %
m_t = zeros(1,N);
m_t(t>0 & t<6)= cos(2*pi*t(t>0 & t<6));
figure;
plot(t,m_t),axis([-2 7]);
title('Signal m(t)');
xlabel('t');
ylabel('x(t)');
grid on
box off

% Step_7.2 Derive an analytical expression for its Fourier transform%
M_f_analytic = 3*sinc(6*(1-f)) + 3*sinc(6*(1+f));

% Step_7.3 Plot analytical expression and Spectrum of sampled signal%
M_f = (fftshift(fft(m_t)))*ts; % dt non-periodic signal
figure;
plot(f,abs(M_f_analytic),'-;analytic expression;','LineWidth',2, f,abs(M_f),'-;Spectrum of sampled signal;','LineWidth',2), axis([-5 5]);
title('|M(f)|');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on
box off

% Step_7.4 Estimate the BW @95% power %
Power_max = max(abs(M_f).^2); %% Or  sum(abs(m).^2)*ts;
Index = find(f == 0);%% Index of zero frequency
for c_index = length(f):-1:length(f)/2
Power = abs(M_f(c_index).^2);
  if(Power >= 0.05*0.5*Power_max)
    BW_m = f(c_index);
    break
  end
end

% Step_8:FDM %
%S1%
fc = 20;
c1_t = cos(2*pi*fc*t);
s1_t = Y1_t.*c1_t;%DSB-SC = Filtered(x_t @ 1Hz).*c1_t
S1_f = fftshift(fft(s1_t))*ts;

%S2%
%% SSB signal
% Step_10 appropriate value for c2(t)%
fc2 = (fc + 2.5 + 1.5);
c2_t = cos(2*pi*fc2*t); % Step 10 %
DSBSC2_t = m_t.*c2_t;% DSB-SC
DSBSC2_f = fftshift(fft(DSBSC2_t))*ts;% frequency spectrum

% Step_9 We are using USB %
BW_channel = 3;
H = zeros(size(f)); %% USB Band Pass Filter
H(f >= fc2 & f <= (fc2 + BW_channel))= 1;%% +ve frequency
H(f <= -fc2 & f >= -(fc2+ BW_channel))= 1;%% +ve frequency
S2_f = DSBSC2_f.*H;%% Filter spectrum
s2_t = ifft(ifftshift(S2_f) / ts);

% Step_11 s_t = s1_t + s2_t %
s_t = s1_t + s2_t;
S_f = fftshift(fft(s_t))*ts;
%s(t) plot
figure;
plot(t,s_t) , axis([-3 7]);
title('s(t)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on
box off
%|S(f)| plot
figure;
plot(f,abs(S_f)) , axis([-28 28]);
title('|S(f)|');
xlabel('Frequency(Hz)');
ylabel('Amplitude');
grid on
box off

%Step_12 coherent demodulator for each channel%
% BPF x(t)%
H1_BP = zeros(size(f));
H1_BP(f > (fc - BW_channel) & f < (fc + BW_channel)) = 1; %% +ve frequency
H1_BP(f < -(fc - BW_channel) & f > -(fc + BW_channel)) = 1; %% -ve frequency
V1_f = S_f .*H1_BP;
v1_t = ifft(ifftshift(V1_f) / ts);

%multiply by the carrier
g1_t = v1_t.*c1_t;
G1_f = fftshift(fft (g1_t))*ts;%ts non-periodic signal

%% LPF (ideal)
H1_LP = abs(f) < BW_channel;
x_t_rec = real(ifft(ifftshift(H1_LP.*G1_f) /ts));

%plot x(t) Received and Transmitted%
figure;
plot(t,x_t_rec./max(x_t_rec))
hold all
plot(t,Y1_t/max(Y1_t)), axis([-5 5]) ;
title('x(t)');
legend('Received message','Transmitted message')
xlabel('Time (s)')
ylabel('Amplitude')
box off

% BPF m(t)%
BW_channel = 3;
H2_BP = zeros(size(f));
H2_BP(f > (fc2) & f < (fc2+ BW_channel)) = 1; %% +ve frequency
H2_BP(f < -(fc2) & f > -(fc2+ BW_channel)) = 1; %% -ve frequency
V2_f = S_f .*H2_BP;
v2_t = ifft(ifftshift(V2_f) / ts);

%multiply by the carrier
g2_t = v2_t.*c2_t;
G2_f = fftshift(fft (g2_t)) *ts;%ts non-periodic signal

%% LPF (ideal)
H2_LP = abs(f) < BW_channel;
m_t_rec = real(ifft(ifftshift(H2_LP.*G2_f) /ts));

%plot m(t) Received and Transmitted%
figure;
plot(t,m_t_rec./max(m_t_rec))
hold all
plot(t,m_t/max(m_t)), axis([-2 8]) ;
title('m(t)');
legend('Received message','Transmitted message')
xlabel('Time (s)')
ylabel('Amplitude')
box off

