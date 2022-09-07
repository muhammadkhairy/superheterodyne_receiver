clear;
clf;
Trange = [0 5 -2 2];
Frange1 = [-100*10^3 100*10^3 0 4000];
Frange2 = [-100*10^3 100*10^3 0 6000];
%% Message
B=6*10^3;
instr ={'select one audio'};
list={'channel 1','channel 2','channel 3'};
result= listdlg('PromptString',instr,'ListString',list,'SelectionMode','single');

[music1, fs] = audioread('sample1.wav');
[music2, fs] = audioread('sample2.wav');
[music3, fs] = audioread('sample3.wav');
info1 = audioinfo('sample1.wav');
info2 = audioinfo('sample2.wav');
info3 = audioinfo('sample3.wav');
durations = [info1.Duration, info2.Duration, info3.Duration];
ts=1/fs;                                              %sampling time              
t=[0 : ts : max(durations)]; t=t(1:end-1);            %time period
L=length(t);
m1 = transpose(music1(:,1)); m1 = [m1, zeros(1, L-length(m1))];
m2 = transpose(music2(:,1)); m2 = [m2, zeros(1, L-length(m2))];
m3 = transpose(music3(:,1)); m3 = [m3, zeros(1, L-length(m3))];
freqs=(-L/2:L/2-1)/(L*ts);

% filtering m to 6KHZ
[bs, as] = butter(6, B/(fs/2));
m1 = 3*filter(bs,as,m1); %Remove the neighboring stations
m2 = 3*filter(bs,as,m2);
m3 = 3*filter(bs,as,m3);

M1=fftshift(fft(m1, L));
M2=fftshift(fft(m2, L));
M3=fftshift(fft(m3, L));

%% Carrier
fc1=40*10^3;                  
c1=cos(2*pi*fc1*t);

fc2=60*10^3;
c2=cos(2*pi*fc2*t);

fc3=80*10^3;
c3=cos(2*pi*fc3*t);

%% Modulation
dsb1 = m1.*(cos(2*pi*fc1*t));
dsb2 = m2.*(cos(2*pi*fc2*t));
dsb3 = m3.*(cos(2*pi*fc3*t));
dsb  = dsb1 + dsb2 + dsb3;

DSB1 = fftshift(fft(dsb1,L));
DSB2 = fftshift(fft(dsb2,L));
DSB3 = fftshift(fft(dsb3,L));
DSB  = fftshift(fft(dsb,L));

%% Demoduation test
dem1 = dsb1 .* (2*cos(2*pi*fc1*t));
dem2 = dsb2 .* (2*cos(2*pi*fc2*t));
dem3 = dsb3 .* (2*cos(2*pi*fc3*t));
% low pass filter
[b,a]=butter(6, B/(fs/2));
dem1=filter(b,a,dem1);
dem2=filter(b,a,dem2);
dem3=filter(b,a,dem3);

DEM1 = fftshift(fft(dem1, L));
DEM2 = fftshift(fft(dem2, L));
DEM3 = fftshift(fft(dem3, L));

%sound(s_dem, fs);

%% plotting
figure(1);
% signal
subplot(3,2,1);
plot(t,m1);
axis(Trange);
title('message signal (t-domain)');

subplot(3,2,2);
fd1 = plot(freqs, abs(M1));
axis(Frange1);
set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('message spectrum (f-domain)');

% modulated
subplot(3,2,3);
plot(t,dsb1);
axis(Trange);
title('modulated signal (t-domain)');

subplot(3,2,4);
fd1 = plot(freqs, abs(DSB1));
axis(Frange1);
set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('modulated signal spectrum (f-domain)');

% demodulated
subplot(3,2,5);
axis(Trange);
plot(t, dem1);
title('De`-modulated with receiver signal (t-domain)');

subplot(3,2,6);
fd1 = plot(freqs, abs(DEM1));
axis(Frange1); 
set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('De-modulated signal spectrum (f-domain)');

figure(2);
% signal
subplot(3,2,1);
plot(t,m2);
axis(Trange);
title('message signal (t-domain)');

subplot(3,2,2);
fd1 = plot(freqs, abs(M2));
axis(Frange1);
set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('message spectrum (f-domain)');

% modulated
subplot(3,2,3);
plot(t,dsb2);
axis(Trange);
title('modulated signal (t-domain)');

subplot(3,2,4);
fd1 = plot(freqs, abs(DSB2));
axis(Frange1);
set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('modulated signal spectrum (f-domain)');

% demodulated
subplot(3,2,5);
axis(Trange);
plot(t, dem2);
title('De`-modulated with receiver signal (t-domain)');

subplot(3,2,6);
fd1 = plot(freqs, abs(DEM2));
axis(Frange1); 
set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('De-modulated signal spectrum (f-domain)');

figure(3);
% signal
subplot(3,2,1);
plot(t,m3);
axis(Trange);
title('message signal (t-domain)');

subplot(3,2,2);
fd1 = plot(freqs, abs(M3));
axis(Frange1);
set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('message spectrum (f-domain)');

% modulated
subplot(3,2,3);
plot(t,dsb3);
axis(Trange);
title('modulated signal (t-domain)');

subplot(3,2,4);
fd1 = plot(freqs, abs(DSB3));
axis(Frange1);
set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('modulated signal spectrum (f-domain)');

% demodulated
subplot(3,2,5);
axis(Trange);
plot(t, dem3);
title('De`-modulated with receiver signal (t-domain)');
%sound(dem3, fs);

subplot(3,2,6);
fd1 = plot(freqs, abs(DEM3));
axis(Frange1); 
set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('De-modulated signal spectrum (f-domain)');
%% Superheterodyne receiver
% if >= B
% fc >= if + B
f_if = 20*10^3;

if(result == 1)
        fc = fc1;
elseif(result == 2 )  
        fc = fc2; 
elseif(result == 3)   
        fc = fc3;

end
% plotting
figure(4);
% plot message
subplot(5,2,1);
plot(t, dsb);
axis(Trange);
title('received modulated message signal (t-domain)');

subplot(5,2,2);
fd1 = plot(freqs, abs(DSB));
axis(Frange2); 
set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('received modulated message spectrum (f-domain)');

%% RF amplifier & filter
[brf, arf] = butter(6, [(fc-B), (fc+B)]/(fs/2));
rfamp = 2*filter(brf, arf, dsb); %Remove the neighboring stations

subplot(5,2,3);
plot(t, rfamp);
axis(Trange);
title('RF Filtered signal (t-domain)');

RFAMP_s=fftshift(fft(rfamp, L));
subplot(5,2,4);
fd1 = plot(freqs, abs(RFAMP_s));
axis(Frange2);
set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('RF Filtered signal spectrum (f-domain)');

%sound(rff,fs);
%% mixer
f_lo = f_if + fc;   %local oscillator frequency
f_gen=cos(2*pi*f_lo*t);
mixed=rfamp.*f_gen;

subplot(5,2,5);
plot(t, mixed);
axis(Trange);
title('Mixed signal (t-domain)');

MIXED = fftshift(fft(mixed,L));

subplot(5,2,6)
title('Response after mixing signals ');
xlabel('Frequency in Hz');
ylabel('Amplitude');
plot(freqs, abs(MIXED));
axis(Frange2);
title('Mixed signal spectrum (f-domain)');

%% IF bandpass filter
lif_cutoff_cof= (f_if-B)/(fs/2); %removing lower adjacent station
uif_cutoff_cof= (f_if+B)/(fs/2); %removing upper adjacent station
[bif, aif] = butter(6,[lif_cutoff_cof, uif_cutoff_cof]);
iff =2* filter(bif,aif,mixed); %Remove the neighboring stations

subplot(5,2,7);
axis(Trange);
plot(t, iff);
axis(Trange);
title('IF Filtered signal (t-domain)');

%cheking in frequency domain
IFF = fftshift(fft(iff, L)); %Freq Response of Filtered Stations

subplot(5,2,8)
plot(freqs, abs(IFF));
axis(Frange2);
title('IF Filter Response (f-domain) ')
xlabel('Frequency (Hz)')
ylabel('Magnitude')

%% receiver demodulation
recdem = iff .* (2*cos(2*pi*(f_if)*t));
%%s_recdem = s_recdem .* (2*cos(2*pi*(fc)*t));
% low pass filter
cutoff_freq = 8000;
cutoff_coefficient = cutoff_freq/(fs/2);
[b,a]=butter(5, cutoff_coefficient);
recdem=filter(b,a,recdem);

RECDEM = fftshift(fft(recdem, L));

sound(recdem, fs);

subplot(5,2,9);
plot(t, recdem);
axis(Trange);
title('Demodulated signal with receiver signal (t-domain)');

%cheking in frequency domain, verifying fft
subplot(5,2,10);
plot(freqs, abs(RECDEM));
axis(Frange2);
title('Demodulated signal with receiver signal (f-domain)');
xlabel('Frequency (Hz;)');
ylabel('Magnitude');
