

%data = mod(randn(2,100000),2*pi); % simulate data for testing

F = 1:150;
nElecs = size(data,1);
CohMat = nan(nElecs,nElecs,length(F));
fs=500;
WINDOW = 1.0 * fs;
NOVERLAP = round(0.75*fs);
%rat = tio';


for e_i = 1:nElecs
    [thetaPhs(e_i,:), thetaAmp(e_i,:)] = extractThetaPhase(data(e_i,:),fs,'hilbert',[5 9]);
end

% 
for e_i = 1:nElecs
    e_i
    for e_j = e_i:nElecs
%         [tmp, F] = mscohere(rat(e_i,:), rat(e_j,:),WINDOW,NOVERLAP,F,fs);
%          CohMat(e_j,e_i,:) = tmp;
%         cohere = (mean(CohMat(:,:,5:9),3));
        
        [plv_th(e_j,e_i), plv_rb(e_j,e_i)] = circmean(circDiff2([thetaPhs(e_j,:); thetaPhs(e_i,:)]),2);
    end
end
figure; imagesc(mean(CohMat(:,:,5:9),3),[0 1]); colorbar; title('theta');
figure; imagesc(mean(CohMat(:,:,30:100),3),[0 1]); colorbar; title('gamma');

