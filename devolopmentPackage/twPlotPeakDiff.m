function [] = twPlotPeakDiff(D)

%need to filter and run through hilbert... 

for iD =  1:size(dataDS,1)
        tInfo.theta_filt(iD,:) = filtfilt(tInfo.btheta,tInfo.atheta,dataDS(iD,:)); %filter the data
        tInfo.theta_phase(iD,:) = atan2(imag(hilbert(tInfo.theta_filt(iD,:))), tInfo.theta_filt(iD,:));
        tInfo.theta_amp(iD,:) = abs(hilbert(tInfo.theta_filt(iD,:)));
end
    
chans = size(D,1);


  for iC1 = 1:length(chans)
    for iC2 = 1:length(chans)
      dataShift(iC1,iC2) = circDiff([D(iC1,:)', D(iC2,:)'],2,'rad'); 
      % changed orientation of data
    end
  end
    

end