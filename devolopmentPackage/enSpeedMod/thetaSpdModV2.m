function [thetaSpdMod_A,thetaSpdMod_F,fullResults] = thetaSpdModV2(self,fRange,useRot)
% function [thetaSpdMod_A,thetaSpdMod_F,fullResults] = thetaSpdModV2(self,fRange,useRot)
% 
% 130906 - eln
% V2 cleaned up some unused code from previous version and
% implemented greater analysis of acceleration as part of core analysis.


  if isempty(self.active_lfp), error('Please set active_lfp'); end

  useUnbinnedData = 0;

  % compute velocity data
  b_velHD = computeSpeed(self,1,useRot);
  
  % compute acceleration data
  b_accHD = diff(b_velHD).*self.fs_video;
  b_accHD = interp1(1:length(b_accHD),b_accHD,linspace(1,length(b_accHD),length(b_accHD)+1));
  

  % quit now if speed data sucks!
  if sum(b_velHD>100)/length(b_velHD)>0.1, % more than 10% of the time he ran more than 100cm/sec (WHOA!)
    [thetaSpdMod_A,thetaSpdMod_F,fullResults] = deal(nan);
    return
  end

  % compute theta power
  if isempty(self.b_lfp(self.active_lfp).signal), self = self.LoadLFP(self.active_lfp,'downsample',250); end
  signal = self.b_lfp(self.active_lfp).signal; fs = self.b_lfp(self.active_lfp).fs;
  %[signal,fs] = scaledEEG(self,0,0);
  fRangeFull = fRange(1)-2:0.125:fRange(end);
  [~,thetaamps] = multiphasevec(fRangeFull,signal,fs,8);
  [thetaAmpl, thetaFreq] = max(thetaamps);
  thetaFreq = fRangeFull(thetaFreq);
  thetaAmpl = log10(thetaAmpl); % use log power to reduce 1/f background
  thetaSig  = buttfilt(signal,[fRange(1) fRange(end)],fs,'bandpass',3);
  
  % match theta fs to vel fs and apply epochs
  b_velHD = interp1(self.b_ts,b_velHD',self.ts)';
  thetaAmpl = interp1(self.b_lfp(self.active_lfp).ts,thetaAmpl',self.ts)';
  thetaFreq = interp1(self.b_lfp(self.active_lfp).ts,thetaFreq',self.ts)';

  % drop nan values
  b_velHD(thetaFreq<fRange(1)) = [];
  thetaAmpl(thetaFreq<fRange(1)) = [];
  thetaFreq(thetaFreq<fRange(1)) = [];

  % quit now if eeg data sucks!
  if length(b_velHD)/length(self.ts) < 0.2
    [thetaSpdMod_A,thetaSpdMod_F,fullResults] = deal(nan);
    return
  end

  % bit more clean up
  velHD = b_velHD(1:length(thetaAmpl));
  accHD = b_accHD(1:length(thetaAmpl));
  thetaAmpl(isinf(thetaAmpl))=nan;

  % define sampling windows for vel and acc
  vel_dim = 0:2.5:max(b_velHD);
  acc_dim = min(b_accHD):5:max(b_accHD);

  if max(thetaAmpl)<100, amp_dim = 0:0.1:3; else amp_dim = 0:10:300; end


  [T_byVel, vInds] = histc(velHD,vel_dim);
  thFrq_byVel = nan(1,length(vel_dim));
  thAmp_byVel = nan(1,length(vel_dim));
  thFrq_byVelHist = nan(length(fRangeFull(fRangeFull>=fRange(1))), length(vel_dim));
  thAmp_byVelHist = nan(length(amp_dim), length(vel_dim));
  [T_byAccVel,thFrq_byAccVel,thAmp_byAccVel] = deal(nan(length(acc_dim),length(vel_dim)));  
  for v = 1:length(vel_dim)
    % bin data by running speed bins
    tmpThetaFrq = thetaFreq(vInds==v);
    tmpThetaAmp = thetaAmpl(vInds==v);
    thFrq_byVel(v) = nanmean(tmpThetaFrq);
    thAmp_byVel(v) = nanmean(tmpThetaAmp);
    thFrq_byVelHist(:,v) = histc(tmpThetaFrq,fRangeFull(fRangeFull>=fRange(1)));
    thAmp_byVelHist(:,v) = histc(tmpThetaAmp,amp_dim);    
    % bin data for vel-by-acc conjunction
    tmpAccHD = accHD(vInds==v);
    [T_byAccVel(:,v), aInds] = histc(tmpAccHD,acc_dim);
    for a = 1:length(acc_dim)
      thFrq_byAccVel(a,v) = nanmean(thetaFreq(aInds==a));
      thAmp_byAccVel(a,v) = nanmean(thetaAmpl(aInds==a));
    end
  end

  
  % bin data into acceleration bins
  [T_byAcc, inds] = histc(accHD,acc_dim);
  thFrq_byAcc = nan(1,length(acc_dim));
  thAmp_byAcc = nan(1,length(acc_dim));
  thFrq_byAccHist = nan(length(fRangeFull(fRangeFull>=fRange(1))), length(acc_dim));
  thAmp_byAccHist = nan(length(amp_dim), length(acc_dim));
  for a = 1:length(acc_dim)
    thFrq_byAcc(a) = nanmean(thetaFreq(inds==a));
    thAmp_byAcc(a) = nanmean(thetaAmpl(inds==a));
    thFrq_byAccHist(:,a) = histc(thetaFreq(inds==a),fRangeFull(fRangeFull>=fRange(1)));
    thAmp_byAccHist(:,a) = histc(thetaAmpl(inds==a),amp_dim);
  end

  % bin frq data by theta amp
  [T_byAmp, inds] = histc(thetaAmpl,amp_dim);
  thFrq_byAmp = nan(1,length(T_byAmp));
  thFrq_byAmpHist = nan(length(fRangeFull(fRangeFull>fRange(1))), length(T_byAmp));
  for a = 1:length(T_byAmp)
    thFrq_byAmp(a) = nanmean(thetaFreq(inds==a));
    thFrq_byAmpHist(:,a) = histc(thetaFreq(inds==a),fRangeFull(fRangeFull>fRange(1)));
  end

  % bin amp data by theta frq
  frq_dim = min(fRange):0.1:max(fRange);
  [T_byFrq, inds] = histc(thetaFreq,frq_dim);
  thAmp_byFrq = nan(1,length(T_byFrq));
  for f = 1:length(T_byFrq)
    thAmp_byFrq(f) = nanmean(thetaAmpl(inds==f));
  end

  
  % set T_* variables to be percent time
  T_byAccVel = T_byAccVel ./ sum(sum(T_byAccVel));
  T_byVel = T_byVel ./ sum(T_byVel);
  T_byAcc = T_byAcc ./ sum(T_byAcc);
  T_byAmp = T_byAmp ./ sum(T_byAmp);
  T_byFrq = T_byFrq ./ sum(T_byFrq);

  
  % choose window to compute slope over running speeds
  enoughVel = velHD>=5 & velHD<=30;
  thetaFreq_ = thetaFreq(enoughVel);
  thetaAmpl_ = thetaAmpl(enoughVel);

  % compute frq vs amp & amp vs Frq slopes
  [frqModAmp(1,:),~,~,~,frqModAmp_R(1,:)] = regress(thetaAmpl_(:), [thetaFreq_(:), ones(size(thetaAmpl_(:)))]);
  [ampModFrq(1,:),~,~,~,ampModFrq_R(1,:)] = regress(thetaFreq_(:), [thetaAmpl_(:), ones(size(thetaAmpl_(:)))]);

  % compute vel vs frq and amp slopes
  if useUnbinnedData % perform regression on all data
    velHD_ = velHD(enoughVel);
    [thetaSpdMod_F(1,:),~,~,~,thetaSpdMod_F_R(1,:)] = regress(thetaFreq_(:),[velHD_(:),ones(size(velHD_(:)))]);
    [thetaSpdMod_A(1,:),~,~,~,thetaSpdMod_A_R(1,:)] = regress(thetaAmpl_(:),[velHD_(:),ones(size(velHD_(:)))]);
  else % compute regression on the binned data
    velDimInds = vel_dim>=5 & vel_dim<=30 & T_byVel>0.005;
    [thetaSpdMod_F(1,:),~,~,~,R] = regress(thFrq_byVel(velDimInds)',[vel_dim(velDimInds);ones(size(vel_dim(velDimInds)))]'); if ~all(isnan(R)), thetaSpdMod_F_R(1,:) = R; else thetaSpdMod_F_R(1,1:4) = nan; end
    [thetaSpdMod_A(1,:),~,~,~,R] = regress(thAmp_byVel(velDimInds)',[vel_dim(velDimInds);ones(size(vel_dim(velDimInds)))]'); if ~all(isnan(R)), thetaSpdMod_A_R(1,:) = R; else thetaSpdMod_A_R(1,1:4) = nan; end
  end


  % now acc vs frq and amp compute slopes
  if useUnbinnedData % perform regression on all data
    enoughAcc = accHD>=-30 & AccHD<=30;
    thetaFreq_ = thetaFreq(enoughAcc);
    thetaAmpl_ = thetaAmpl(enoughAcc);
    accHD_ = accHD(enoughAcc);
    [thetaAccMod_F(1,:),~,~,~,thetaAccMod_F_R(1,:)] = regress(thetaFreq_(:),[accHD_(:),ones(size(accHD_(:)))]);
    [thetaAccMod_A(1,:),~,~,~,thetaAccMod_A_R(1,:)] = regress(thetaAmpl_(:),[accHD_(:),ones(size(accHD_(:)))]);
  else % compute regression on the binned data
    accDimInds = acc_dim>=-30 & acc_dim<=30 & T_byAcc>0.005;
    [thetaAccMod_F(1,:),~,~,~,R] = regress(thFrq_byAcc(accDimInds)',[acc_dim(accDimInds);ones(size(acc_dim(accDimInds)))]'); if ~all(isnan(R)), thetaAccMod_F_R(1,:) = R; else thetaAccMod_F_R(1,1:4) = nan; end
    [thetaAccMod_A(1,:),~,~,~,R] = regress(thAmp_byAcc(accDimInds)',[acc_dim(accDimInds);ones(size(acc_dim(accDimInds)))]'); if ~all(isnan(R)), thetaAccMod_A_R(1,:) = R; else thetaAccMod_A_R(1,1:4) = nan; end
  end


  % pack it all up
  fullResults.T_byVel = T_byVel(T_byVel>0.005);
  fullResults.vel = vel_dim(T_byVel>0.005);
  fullResults.thFrq_byVel = thFrq_byVel(T_byVel>0.005);
  fullResults.thAmp_byVel = thAmp_byVel(T_byVel>0.005);
  fullResults.thetaSpdMod_F = thetaSpdMod_F;
  fullResults.thetaSpdMod_F_R = thetaSpdMod_F_R;
  fullResults.thetaSpdMod_A = thetaSpdMod_A;
  fullResults.thetaSpdMod_A_R = thetaSpdMod_A_R;

  fullResults.T_byAcc = T_byAcc(T_byAcc>0.005);
  fullResults.acc = acc_dim(T_byAcc>0.005);
  fullResults.thFrq_byAcc = thFrq_byAcc(T_byAcc>0.005);
  fullResults.thAmp_byAcc = thAmp_byAcc(T_byAcc>0.005);
  fullResults.thetaAccMod_F = thetaAccMod_F;
  fullResults.thetaAccMod_F_R = thetaAccMod_F_R;
  fullResults.thetaAccMod_A = thetaAccMod_A;
  fullResults.thetaAccMod_A_R = thetaAccMod_A_R;

  fullResults.T_byAccVel = T_byAccVel(T_byAcc>0.005,T_byVel>0.005);
  fullResults.thFrq_byAccVel = thFrq_byAccVel(T_byAcc>0.005,T_byVel>0.005);
  fullResults.thAmp_byAccVel = thAmp_byAccVel(T_byAcc>0.005,T_byVel>0.005);
  
  fullResults.meanAmp = nanmean(thetaAmpl);
  fullResults.amp = amp_dim(T_byAmp>0.005);
  fullResults.T_byAmp = T_byAmp(T_byAmp>0.005);
  fullResults.thFrq_byAmp = thFrq_byAmp(T_byAmp>0.005);
  fullResults.ampModFrq = ampModFrq;
  fullResults.ampModFrq_R = ampModFrq_R;

  fullResults.meanFrq = nanmean(thetaFreq);
  fullResults.frq = frq_dim(T_byFrq>0.005);
  fullResults.T_byFrq = T_byFrq(T_byFrq>0.005);
  fullResults.thAmp_byFrq = thAmp_byFrq(T_byFrq>0.005);
  fullResults.frqModAmp = frqModAmp;
  fullResults.frqModAmp_R = frqModAmp_R;

  fullResults.thAmp_byVelHist = thAmp_byVelHist;
  fullResults.thFrq_byVelHist = thFrq_byVelHist;
  fullResults.thAmp_byAccHist = thAmp_byAccHist;
  fullResults.thFrq_byAccHist = thFrq_byAccHist;
  fullResults.thFrq_byAmpHist = thFrq_byAmpHist;

end




function b_velHD = computeSpeed(self,pureSpeed,useRot)
  if pureSpeed
    % prepare speed data
    self = CMBHOME.Utils.MergeEpochs(self);
    [val,ind] = max(self.epoch);
    if val>max(self.b_lfp(self.active_lfp).ts), self.epoch(ind) = max(self.b_lfp(self.active_lfp).ts); end
    b_velHD = self.b_vel.*self.spatial_scale;
  else
    % project movement onto heading direction to extract lateral versus forward movements
    % prepare speed data
    self = CMBHOME.Utils.MergeEpochs(self);
    [val,ind] = max(self.epoch);
    if val>max(self.b_lfp(self.active_lfp).ts), self.epoch(ind) = max(self.b_lfp(self.active_lfp).ts); end
    b_velVect = [diff(self.b_x), diff(self.b_y)]' .* self.spatial_scale .* self.fs_video;
    [b_velTh,b_velR] = cart2pol(b_velVect(1,:),b_velVect(2,:));
    b_hdTh = deg2rad(self.b_headdir(1:end-1))';

    % update hd to align average velocity and head direction to be the
    % same, assuming some fixed offset in head tracker as mounted to the
    % head stage.
    diffs = circDiff([b_hdTh;b_velTh]',2,'rad');
    trackerOffset = circmean(diffs');
    b_hdTh = mod(b_hdTh+trackerOffset+pi,2*pi)-pi;
    b_hdVect = [cos(b_hdTh); sin(b_hdTh)];
    diffs = circDiff([b_hdTh;b_velTh]',2,'rad');

    % compute velocity to be the projection onto the direction the rat is facing
    [b_hd(:,1),b_hd(:,2)] = pol2cart(diffs,b_velR);

    % decide if we should use rotational velocity intstead of forward velocity
    if ~exist('useRot','var') | isempty(useRot)
      useRot = 0;
    end
    if ~useRot,    b_velHD = b_hd(:,1);
    else b_velHD = b_hd(:,2); end
    b_velHD = interp1(1:length(b_velHD),b_velHD,linspace(1,length(b_velHD),length(b_velHD)+1));
  end
end