

F = 1:150;
nElecs = 5;
CohMat = nan(nElecs,nElecs,length(F));
fs=500;
WINDOW = 1.0 * fs;
NOVERLAP = round(0.75*fs);


for e_i = 1:nElecs-1
    e_i
    for e_j = e_i+1:nElecs
        [tmp, F] = mscohere(rio(e_i,:), rio(e_j,:),WINDOW,NOVERLAP,F,fs);
         CohMat(e_j,e_i,:) = tmp;
    end
end
figure; imagesc(mean(CohMat(:,:,6:9),3),[0 1]); colorbar; title('theta');
figure; imagesc(mean(CohMat(:,:,30:100),3),[0 1]); colorbar; title('gamma');

