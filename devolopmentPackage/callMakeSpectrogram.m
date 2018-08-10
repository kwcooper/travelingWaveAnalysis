
figure; plot(sig(1,1:100000))

x = sig(1,1:100000);
s = spectrogram(x);

figure; 
spectrogram(x,'yaxis')

t = 0:0.001:2;
figure; spectrogram(x,256,50,[],1e3,'yaxis'); ylim([0 100]);

figure; pspectrum(x,500,'spectrogram')

figure;
z = randn(100,100);
t = 1:100;
x = 1:100;
surf(t,x,abs(z),'EdgeColor','none');
axis xy; axis tight; colormap(jet); view(0,90);


[px f] = welch(x, 500, ones(1,500), 250, 100, 0);
figure; imagesc(10*log10(abs(px)))




