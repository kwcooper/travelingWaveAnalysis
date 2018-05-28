

function [thetaSpdMod_R] = thetaSpdAnalysis(root, chans)
% which of the channels should we analyse? Why not all of them! muh hahahah

metaData = root.user_def.metaData;
plt = 1;

% for i = 1:chans
% root.active_lfp = i;
% [thetaSpdMod_amp, thetaSpdMod_frq, thetaSpdMod_R] = thetaSpdMod(root,[6 10],0);
% end

root.active_lfp = 1;
[thetaSpdMod_amp, thetaSpdMod_frq, thetaSpdMod_R] = thetaSpdMod(root,[6 10],0);

t = thetaSpdMod_R;


%%
if plt
  figure; imagesc(t.vel,t.frq,t.thFrq_byVelHist); title(['Freq by velocity: ', metaData.Rat]);
  
  plotName = [metaData.Recording '_' metaData.Rat];
  if ~exist(fullfile(metaData.savePath, 'thetaSpdHist')); mkdir(fullfile(metaData.savePath, 'thetaSpdHist')); end
  printFigure(gcf, [fullfile(metaData.savePath, 'thetaSpdHist',[plotName,'.',metaData.fig_type])],'imgType',metaData.fig_type);
  fprintf('Saved figure (thetaSpdHist)\n');
  
  
  figure; plot(t.vel,t.thAmp_byVel); title(['theta Amp by velocity: ', metaData.Rat]);
  
  plotName = [metaData.Recording '_' metaData.Rat];
  if ~exist(fullfile(metaData.savePath, 'thetaAmp')); mkdir(fullfile(metaData.savePath, 'thetaAmp')); end
  printFigure(gcf, [fullfile(metaData.savePath, 'thetaAmp',[plotName,'.',metaData.fig_type])],'imgType',metaData.fig_type);
  fprintf('Saved figure (thetaAmp)\n');
  
  
  figure; plot(t.vel,t.thFrq_byVel); title(['theta Freq by velocity: ', metaData.Rat]);
  
  plotName = [metaData.Recording '_' metaData.Rat];
  if ~exist(fullfile(metaData.savePath, 'thetaSpd')); mkdir(fullfile(metaData.savePath, 'thetaSpd')); end
  printFigure(gcf, [fullfile(metaData.savePath, 'thetaSpd',[plotName,'.',metaData.fig_type])],'imgType',metaData.fig_type);
  fprintf('Saved figure (thetaSpd)\n');
end
end

