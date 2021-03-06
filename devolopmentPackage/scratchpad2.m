

thetaPhase = root.user_def.theta_phs;

% plot each of the theta phases
figure; suptitle('thetaPhase')
for i = 1:size(thetaPhase,1)
  subplot(size(thetaPhase,1),1,i)
  plot(thetaPhase(i,1:1000), '.');
end

%Finds slope offset (from the first channel)
  for iC1 = 1:size(thetaPhase(1,:),1)
    for iC2 = 1:size(thetaPhase,1)
      dataShift(iC1, iC2, :) = circDiff([thetaPhase(iC1,:)', thetaPhase(iC2,:)'],2,'rad');
      meanDataShift(iC1,iC2) = circmean(circDiff([thetaPhase(iC1,:)', thetaPhase(iC2,:)'],2,'rad'),2); 
    end
  end
 
% grab the first channel diffs into a 2d matrix
s1 = squeeze(dataShift(1,:,:));
w = 10000; figure; plot(s1(:,1:w)'); title('Theta Difference');

ms1 = circmean(s1,1);
w = 10000; figure; plot(ms1(:,1:w)'); title('Mean Theta Difference');


dtp = thetaPhase(3,:) - thetaPhase(4,:);
figure; plot(dtp(:,1:w)');
cdtp = circDiff([thetaPhase(3,:)', thetaPhase(4,:)'],2,'rad');
figure; plot(root.lfp.ts(1:w,:), cdtp(:,1:w)');

q = circDiff([thetaPhase(iC1,:)', thetaPhase(2,:)'],2,'rad');

% 1. tracking with multiple rats
% 2. double check root objects
% 3. fix/add figure saving 
% 4. confirm references/tracking
% 5. build slope vs velocity function
% 6. check timestamps
% 7. check slope gathering functions
% 8. check 



%%
span = linspace(1,10,100);
e = sin(span + 4/pi);
d = sin(span); 
f = sin(span + 2/pi);
figure; hold on;
plot(f); plot(e); plot(d);

emd = d-e; figure; plot(emd);
cdemd = circDiff([e',d'], 2, 'rad'); figure; plot(cdemd);

waveMat = [d; f; e];
figure; 
for i = 1:size(waveMat,1)
  subplot(size(waveMat,1),1,i)
  plot(waveMat(i,:), '.');
end

for i = 1:size(waveMat,1)-1
  diffMat(i,:) = circDiff([waveMat(1,:)',waveMat(i + 1,:)'], 2, 'rad');
end

figure; 
for i = 1:size(diffMat,1)
  subplot(size(diffMat,1),1,i)
  plot(diffMat(i,:), '.');
end

figure; plot(circmean(diffMat, 1)')
figure; plot(mean(diffMat, 1))

%%

 we suffer
 we suffer because we crave
 there exists a cessessation of this suffering
 the eightfold path is the way to the cessation of suffering

 right effort
 right view
 right intention
 right livleyhood
 right speech
 right concentration 
 right action
 right mindfulness
 
 gradual training
 
 desire; ill will; sloth; regret; doubt
 