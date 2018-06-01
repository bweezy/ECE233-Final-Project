%%%%% PROJECT %%%%%%%%%%
clear all; close all; clc;
dbstop if error


%%%%Parameters %%%%%
Nb = 1000;  %Number of bits
bw = 1e9;   %Bandwidth
Fc = 3.5e9; %Carrier Frequency

Nf = 5;    %Number of frames per signal
Nc = 5;    %Number of chips per frame


%%%Derived Parameters%%%%%%%%
f_h = Fc + .5*bw; %High end of bandwidth
f_l = Fc - .5*bw; %Low end of bandwidth

Fs_tx = 100e9; %Sampling rate for tx
Fs_rx = 2*bw;  %Sampling rate for rx subsampling

Tc = 1/bw; %Chip Time
Tf = Tc*Nc; %Frame Time
Ts = Tf*Nf; %Symbol Time
% Ts = 1/bw; %Symbol time
% Tf = Ts/Nf; %Frame time
% Tc = Tf/Nc; %Chip Duration

frame_s = int32(Tf*Fs_tx); %Samples per frame
chip_s = int32(Tc*Fs_tx); %Samples per chip
symbol_s = int32(Ts*Fs_tx); % Samples per symbol


cj = randi([0, Nc-1],1,Nf);  %pseudo-random integer sequence {0, Nc-1}
bits = randi([0,1], 1,Nb);       %Bits (BPSK) {-1, +1}
bits(bits==0) = -1;



%w_tx(t) is the unit-energy transmit waveform
w_tx = gauspuls((-chip_s/2:1:chip_s/2-1),.4,1);
                       
% tc = gmonopuls('cutoff',Fc);      % width of each pulse
% t1  = -2*tc : 1/Fs_tx : 2*tc;
% pulse = gmonopuls(t1,Fc);
t=-5e-10:1/Fs_tx:5e-10;
a= Tc/2.5;
pulse=(1-(4*pi.*(t.^2))/a^2) .* exp(-2*pi.*(t.^2)/a^2) / sqrt(1);

t  = 0 : 1/fs : Nb*Ts*1e-9;
s_tx = pulstran(t, 

s_tx = zeros(1,Nb*symbol_s);

t = (1:double(Nb*symbol_s))/Fs_tx;



for b = 0:Nb-1
   for j = 0:Nf-1
      for c = 0:Nc-1
          if cj(j+1) == c
             s_tx(b*symbol_s+j*frame_s+c*chip_s+1:b*symbol_s+j*frame_s+(c+1)*chip_s) = bits(b+1)*w_tx; 
          else
             s_tx(b*symbol_s+j*frame_s+c*chip_s+1:b*symbol_s+j*frame_s+(c+1)*chip_s) = zeros(1,chip_s);
          end
      end
   end
end

s_tx = s_tx .* cos(2*pi*Fc*t);
Y = fft(s_tx);
L = length(s_tx);

Y = fft(s_tx);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs_tx*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')



