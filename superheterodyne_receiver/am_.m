clear;
clf;
Trange = [0 5 -2 2];
Frange1 = [-100*10^3 100*10^3 -2000 2000];
Frange2 = [-100*10^3 100*10^3 -4000 4000];
%% Message
B=8*10^3;

[music, fs] = audioread('sample1.wav');
info = audioinfo('sample1.wav');
ts=1/fs;                                             %sampling time              
t=[0 : ts : info.Duration]; t=t(1:end-1);            %time period
m = transpose(music(:,1));
L=length(t);
freqs=(-L/2:L/2-1)/(L*ts);

% filtering m to 8KHZ
[bs, as] = butter(6, B/(fs/2));
m = 3*filter(bs,as,m); %Remove the neighboring stations

figure(1);
subplot(3,2,1);
plot(t,m);
axis(Trange);
title('message signal (t-domain)');


M=fftshift(fft(m, L));
subplot(3,2,2);
fd1 = plot(freqs, abs(M));
axis(Frange1);
set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('message spectrum (f-domain)');


%% Carrier
fc=36*10^3;                  %carrier frequency
c=cos(2*pi*fc*t);            %carrier function


% subplot(8,2,3);
% plot(t,c);
% 
% title('carrier signal');
% 
% 
% C_s=fftshift(fft(c, L));
% Frange=[-20*10^6 20*10^6 -1*10^3 1*10^3];
% subplot(8,2,4);
% fd1 = plot(freqs, abs(C_s));
% axis(Frange); set(fd1,'Linewidth',1.5);
% xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
% title('carrier spectrum (f-domain)');



%% Modulation
dsb=m.*(cos(2*pi*fc*t));

subplot(3,2,3);

plot(t,dsb);
axis(Trange);
title('modulated signal (t-domain)');

DSB=fftshift(fft(dsb,L));

subplot(3,2,4);
fd1 = plot(freqs, abs(DSB));
axis(Frange1);
set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('modulated signal spectrum (f-domain)');

%% Demoduation
dem = dsb .* (2*cos(2*pi*fc*t));
% low pass filter
cutoff_freq = 8000;
cutoff_coefficient = cutoff_freq/(fs/2);
[b,a]=butter(5, cutoff_coefficient);
dem=filter(b,a,dem);
S_dem = fftshift(fft(dem, L));


subplot(3,2,5);
axis(Trange);
plot(dem);
title('De`-modulated with receiver signal (t-domain)');

DEM=fftshift(fft(dem, L));
subplot(3,2,6);
fd1 = plot(freqs, abs(DEM));
axis(Frange1); 
set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('De-modulated signal spectrum (f-domain)');

%sound(s_dem, fs);

%% Superheterodyne receiver
% if >= B
% fc >= if + B
f_if = 18*10^3;
fcmax = 74*10^3;



% plotting
figure(2);
% plot message
subplot(6,2,1);
plot(t, m);
axis(Trange);
title('message signal (t-domain)');

subplot(6,2,2);
fd1 = plot(freqs, abs(M));
axis(Frange2); 
set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('message spectrum (f-domain)');

% plot modulated
subplot(6,2,3);
plot(t, dsb);
axis(Trange);
title('modulated signal (t-domain)');

subplot(6,2,4);
fd1 = plot(freqs, abs(DSB));
axis(Frange2); 
set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('modulated signal spectrum (f-domain)');

%% RF filter
rfamp = 3*dsb; %Remove the neighboring stations

subplot(6,2,5);
plot(t, rfamp);
axis(Trange);
title('RF Filtered signal (t-domain)');

RFAMP_s=fftshift(fft(rfamp, L));
subplot(6,2,6);
fd1 = plot(freqs, abs(RFAMP_s));
axis(Frange2);
set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('RF Filtered signal spectrum (f-domain)');

%sound(rff,fs);
%% mixer
f_lo = f_if + fc;   %local oscillator frequency
f_gen=cos(2*pi*f_lo*t);
mixed=2.*rfamp.*f_gen;

subplot(6,2,7);
plot(t, mixed);
axis(Trange);
title('Mixed signal (t-domain)');

MIXED_s = fftshift(fft(mixed,L));

subplot(6,2,8)
title('Response after mixing signals ');
xlabel('Frequency in Hz');
ylabel('Amplitude');
plot(freqs, abs(MIXED_s));
axis(Frange2);
title('Mixed signal spectrum (f-domain)');

%% IF bandpass filter
lif_cutoff_cof= (f_if-B)/(fs/2); %removing lower adjacent station
uif_cutoff_cof= (f_if+B)/(fs/2); %removing upper adjacent station
[bif, aif] = butter(6,[lif_cutoff_cof, uif_cutoff_cof]);
iff = 2*filter(bif,aif,mixed); %Remove the neighboring stations

subplot(6,2,9);
axis(Trange);
plot(t, iff);
title('IF Filtered signal (t-domain)');

%cheking in frequency domain
IFF_s = fftshift(fft(iff, L)); %Freq Response of Filtered Stations

subplot(6,2,10)
plot(freqs, abs(IFF_s));
title('IF Filter Response (f-domain) ')
xlabel('Frequency (Hz)')
ylabel('Magnitude')

%% receiver demodulation
s_recdem = iff .* (2*cos(2*pi*(f_if)*t));
%%s_recdem = s_recdem .* (2*cos(2*pi*(fc)*t));
% low pass filter
cutoff_freq = 8000;
cutoff_coefficient = cutoff_freq/(fs/2);
[b,a]=butter(5, cutoff_coefficient);
s_recdem=filter(b,a,s_recdem);

S_recdem = fftshift(fft(s_recdem, L));

sound(s_recdem, fs);

subplot(6,2,11);
plot(t, s_recdem);
axis(Trange);
title('Demodulated signal with receiver signal (t-domain)');

%cheking in frequency domain, verifying fft
subplot(6,2,12);
plot(freqs, abs(S_recdem));
title('Demodulated signal with receiver signal (f-domain)');
xlabel('Frequency (Hz;)');
ylabel('Magnitude');
