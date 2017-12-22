phases = angle(theta_phase);
band = [6,10]
[cycles,~] = parseThetaCycles(phases,root.fs_video,band); %eln add 170216
inds = find(cycles);