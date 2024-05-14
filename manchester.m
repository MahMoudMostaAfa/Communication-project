clc;
clear all;
close all;
%%%% MANCHESTER Pulses in Time Domain %%%%
m = randi([0 1],1,64);
%%64 bits random bits generator is just a random integer generator of IMIN = 0 , IMAX = 1, 64 row elements and 1 column
T=1000; %total sampling time of signal
df=1/T;
fs = 1000 ; %sampling rate
ts = 1/fs ; %sampling time
N = ceil(T/ts); #number of samples
t = 0 : ts : ((N-1)*ts); #time vector
n = length(m);
x = [];
y = [];
max = 1;
min = -1;
for i=1:n
    x=[x i-1 i-1+0.5 i-1+0.5 i];  %timing vector % x = [0 0.5 0.5 1]
    if(m(i) == 1)   % y == Symbol ‾|_
        y=[y max max min min]; %Data vector % y = [1.2 1.2 -1.2 -1.2]
    else            %(m(i) == 0) %y == Symbol _|⁻
        y=[y min min max max];
    end
end

figure(1);
title('Time Domain of Manchester Pulse')
plot(x,y),
axis([0,n,-2,2]);
grid on
box off
ylabel('Signal(V)');
xlabel('Time(s)');
title('Manchester');



%MANCHESTER Pulses in Frequency Domain

y_zeroFilling = [y, zeros(1, N - length(y))];  %filling the rest of Y with zero after the real signal to improve the resolution
Y_polar = (fftshift(fft(y_zeroFilling)) ) *ts ; %getting the Fourier transform of y
if(rem(N,2)==0)
  f = - (0.5*fs) : df : (0.5*fs-df) ;
else
  f = - (0.5*fs-0.5*df) : df : (0.5*fs-0.5*df);
end

figure(2);
plot(f, abs(Y_polar));
grid on
box off
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Domain Spectrum of Manchester Pulse');

%%%% End of Code %%%%

