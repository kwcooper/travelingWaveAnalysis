function twCombineFigs(figData)

figure;

%Raw waves Fig
subplot(2,2,1);
plot(figData.rawWaves.t,figData.rawWaves.lfpO,'k'); axis ij;
xlabel('Time') %!! Is this correct? or should it be phase?
ylabel('Channel') % !! What about the axis though...
title('Raw Waves')
ltr = text(0.02,0.98,'A','Units', 'Normalized', 'VerticalAlignment', 'Top');
s = ltr.FontSize; ltr.FontSize = 12;
grid on


subplot(2,2,2);



%Average Waves Plot
subplot(2,2,3);
plot(figData.CTA.t,figData.CTA.lfp_,'k'); axis ij;
xlabel('Time') %!! Is this correct? or should it be phase?
ylabel('Channel') % !! What about the axis though...
title('Averaged Waves')
ltr = text(0.02,0.98,'C','Units', 'Normalized', 'VerticalAlignment', 'Top');
s = ltr.FontSize; ltr.FontSize = 12;
grid on;

subplot(2,2,4);

end