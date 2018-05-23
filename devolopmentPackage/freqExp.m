% (D, fs)
% play with other frequency bands


dSamp = D(1, 1:10000);
fs = 500;
figure; plot(dSamp);

%             dl     th       bt       lg       mg
freqList = {[1 4], [6 12], [10 25], [25 60], [60 120]};
freqNames = {'dl', 'th', 'bt', 'lg', 'mg'};
theta = [freqList{2}];

figure; suptitle('Freq Bands')
for i = 1:size(freqList,2)
  subplot(size(freqList,2),1,i)
  plot(buttfilt(dSamp,freqList{i},fs,'bandpass',3))
  title(freqNames{i})
end

figure; suptitle('Hilbert')
for i = 1:size(freqList,2)
  subplot(size(freqList,2),1,i)
  plot(abs(hilbert(buttfilt(dSamp,freqList{i},fs,'bandpass',3))))
  title(freqNames{i})
end






% in bins
f = 0:length(dSamp)-1;
figure; plot(f,abs(fft(dSamp)))


data = buttfilt(dSamp(1:1000),[6 20],fs,'bandpass',3);
data = dSamp(1:1000);
% in hz
% magnitude
fHz = (0:(length(data)-1)) * fs / length(data);
figure; plot(fHz, abs(fft(data)))
xlim([0 100])

% power
Y = abs(fft(data));
fHz = (0:(length(data)-1)) * fs / length(data);
yPow = 10*log10(Y(1:ceil(length(data)/2)));
figure; plot(fHz(1:ceil(length(data)/2)), yPow)

%%


% grab mag for all samples
 data = dSamp;
 maxFreq = 100;
 wSize = 100;
 w = floor(size(data,2) / wSize);
 freqMag = nan(maxFreq, w); 
 id = 1;
 for i = 1:wSize:size(data,2)
   Y = abs(fft(data(i:i+wSize-1)))';
   freqMag(:,id) = Y;
   id = id + 1;
 end
 figure; plot(fHz,freqMag)
 figure; imagesc(freqMag/max(freqMag))

 
 
 
 % grab pow for all samples
 data = D(1,1:700000);
 %data = buttfilt(D(1,1:700000),[25 60],fs,'bandpass',3);

 wSize = fs;
 w = floor(size(data,2) / wSize);
 freqPow = nan(100, w); 
 id = 1;
 for i = 1:wSize:size(data,2)
   yPow = 10*log10(abs(fft(data(i:i+wSize-1))))';
   freqPow(:,id) = yPow(1:100);
   id = id + 1;
 end
 fHz = (0:(w-1)) * fs / w;
 %figure; plot(fHz(1:100),freqPow)
 figure; imagesc(freqPow)
 title('PSD: Tio'); xlabel('Time (s)'); ylabel('Hz');

%% graveyard

% m = length(dSamp2); 
% n = pow2(nextpow2(m));
% y = fft(dSamp2,n); 
% f = (0:n-1)*(fs/n)/10; % frequency vector
% power = abs(y).^2/n;   % power spectrum      
% figure; plot(f(1:floor(n/2)),power(1:floor(n/2)))
% xlabel('Frequency')
% ylabel('Power')


% Y = fft(dSamp);
% n = length(dSamp);
% yShift = fftshift(Y);
% fshift = (-n/2:n/2-1)*(50/n);
% power = abs(yShift).^2/n; 
% figure; plot(fshift,power)
% title('Signal Power')
