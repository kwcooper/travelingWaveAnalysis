function [] = plotPSD(data,fs,tText)
%  plots the Power spectral Density over windows of data
%  (FFT over windows of data)
%
%  wSize - currently set to the sample rate
%  
%  K - 180521

% TD: this only computes one channel, how does it change over the whole
% hippocampus?


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
 title(tText); xlabel('Time (s)'); ylabel('Hz');

end