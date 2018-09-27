% requires ratsComp (NOV FAM)

% Jeewajee 2009
% EEG data were filtered to admit only portions of EEG where the rat 
% spent ? 0.5 seconds at speeds above 5 cm/s.
% 
% constraint: that the median speed in each trial was constant and equal 
% to the median speed of the rat across all trials in the experiment k. 
% Thus, speed limits s1 > 5 cm/s and s2 < ? were chosen such that the length 
% of ordered data was symmetric (and maximal) about k. 
%
% Power spectra were calculated by finding the fast Fourier transform 
% of the concatenated filtered data, where the square-modulus of each 
% Fourier frequency coefficient represents the signal power at that frequency. 
% Finally the power spectrum was smoothed using a Gaussian kernel
% with standard deviation 0.375Hz.
