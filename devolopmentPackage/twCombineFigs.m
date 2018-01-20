function twCombineFigs(figData, plt)

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

%peak offset
subplot(2,2,2);
plot(figData.pd.x,figData.pd.dataShift,'o');
hold on;
plot(figData.pd.x,figData.pd.y);
title([figData.ratInfo.name ' peak offset | Slope: ' num2str(figData.pd.B(1))]);


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



if plt
  plotName = [figData.ratInfo.recording '_' figData.ratInfo.name];
  printFigure(gcf, [fullfile(figData.savePath, 'subPlots',[plotName,'.',figData.fig_type])],'imgType',figData.fig_type);
  fprintf('Saved subplot\n');
end
end