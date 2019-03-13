
function [plv_th, plv_rb] = phaseLocking(data,fs)

%data = mod(randn(2,100000),2*pi); % simulate data for testing
%data = sig;
%fs=500;
%F = 1:150;
nElecs = size(data,1);

% Grab theta phase for each channel
for e_i = 1:nElecs
    [thetaPhs(e_i,:), thetaAmp(e_i,:)] = extractThetaPhase(data(e_i,:),fs,'hilbert',[5 9]);
end

% compute the difference between phases
% keeping the circular mean difference and mean resultant length 
for e_i = 1:nElecs
    for e_j = e_i:nElecs
       [plv_th(e_j,e_i), plv_rb(e_j,e_i)] = circmean(circDiff2([thetaPhs(e_j,:); thetaPhs(e_i,:)]),2);
    end
end

end


