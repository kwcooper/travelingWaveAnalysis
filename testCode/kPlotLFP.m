function h = kPlotLFP(lfp,inds,fs,scale,h)
%function h = plotLFP(lfp,fs,scale,h,spacebarScrolling,useImg)
% 
if ~exist('h','var')
  h = figure;
end

if ~exist('scale','var')
  scale = 2.5;
end

if ~exist('fs','var')
  fs = 1;
end

if ~exist('spacebarScrolling','var')
  spacebarScrolling = 1;
end

if ~exist('useImg','var')
  useImg = 0;
end

[nElecs,tPts,nSets] = size(lfp);

% prepare for multiset plotting if needed
nR = floor(sqrt(nSets));
nC = ceil(nSets/nR);

t = (1:tPts)/fs;
disp(scale)
lfp_ = lfp / (-1 * scale * rms(lfp(:))); % -1 is used to correct for sign flip that comes with use of axis ij below
if ~useImg
  offsets = repmat([1:nElecs]',1,tPts,nSets);
  lfp_ = lfp_ + offsets;
end
%csd = diff(lfp,2);

figure(h); 
%imagesc(t,[2:nElecs-1],csd); hold on;
if ~spacebarScrolling
  plot(t,lfp_);
  % keyboard
  
  xlim([0 1]);
  ylim([0 31])
end


% else
%     fprintf('Press spacebar to continue (ctrl+c to esc)\n');
%     for i = 1:max(1,floor(size(lfp_,2)/fs))
        for s = 1:nSets
            subplot(nR,nC,s);
            if useImg
                imagesc(t(inds),1:size(lfp_,1),lfp_(:,inds,s),[-2 2]);axis ij
            else
                plot(t(inds),lfp_(:,inds,s));axis ij
            end
            ylim([-1 nElecs+1]); grid on;
        end
%         pause
%     end
end

