fs = 500;
thetaRange = [6 10];


lfp = hpcD(1,:);

% grap theta Amp for whole recording
[~,thetaAmp,~] = extractThetaPhase(lfp,fs,'hilbert',thetaRange);
thetaAmp = thetaAmp';
figure; plot(thetaAmp(:,1:100000));

% Theta delta ratio (use greater than 2)  (Sirota et al., 2008)
th_Amp = abs(hilbert(buttfilt(lfpSnipRaw,[6 12],fs,'bandpass',4)));
dl_Amp = abs(hilbert(buttfilt(lfpSnipRaw,[2 4],fs,'bandpass',4)));
TDRatio = th_Amp ./ dl_Amp;

[~,incInds2] = CMBHOME.Utils.OverThresholdDetect(TDRatio,filtParams(1),ceil(fs/4),ceil(fs/2));
excInds = ~incInds2;



% grab a snip of the data
os = 55000;
ws = 500; 
lfpSnipRaw = lfp(1,1+os:ws+os);
figure; plot(lfpSnipRaw);

lfpSnipTheta = buttfilt(lfpSnipRaw',thetaRange,fs,'bandpass',3); % filter for theta
hold on; plot(lfpSnipTheta)

lfpSnipBroad = buttfilt(lfpSnipRaw',[1 80],fs,'bandpass',3); % belluscio 2012
hold on; plot(lfpSnipBroad)
legend({'Raw LFP','6-10Hz filter','1 80Hz filter'})

% find theta phase (peaks?) (Le Van Quyen et al., 2001)
thetaPhs = angle(hilbert(buttfilt(lfpSnipRaw,thetaRange,fs,'bandpass',3))) + pi;
figure; plot(lfpSnipBroad/max(lfpSnipBroad)); hold on; plot(thetaPhs);

%%
% let's play with wavelets
% sin * gauusian = Morlet wavelet

% parameters
sampleRate = 500; 
waveletFreq = 6.5; % wavelet frequency in Hz
%t = linspace(0,1,100);
t = -2:1/sampleRate:2;

sinWave = sin(2*pi*freq*t);
%sinComplex = exp( 1i*2*pi*freq.*time );

s = 7 / (2*pi*waveletFreq); 
gauss = exp(-t^2 / (2 * s^2));

wavelet = gauss .* sinwave; 


%%
% find where there is high theta power


%%
x = linspace(0,1,1000);

base = 4*cos(2*pi*x);

Pos = [1 2 3 5 7 8]/10;
Hgt = [3 7 5 5 4 5];
Wdt = [1 3 3 4 2 3]/100;

for n = 1:length(Pos)
    Gauss(n,:) =  Hgt(n)*exp(-((x - Pos(n))/Wdt(n)).^2);
end

PeakSig = sum(Gauss)+base;

plot(x,Gauss,'--',x,PeakSig,x,base)

figure; findpeaks(PeakSig,x,'MinPeakProminence',4,'Annotate','extents')


%findpeaks(select,Fs,'MinPeakDistance',0.005)
%findpeaks(select,Fs,'MinPeakHeight',1)














