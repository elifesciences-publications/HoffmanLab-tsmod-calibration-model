%% Script to make rc to E conversion lookup table
% Lookup table describes the relationship between (r/ro)^6 and
% FRET efficiency in the static random isotropic averaging regime
clear;
for ro = 4.5:0.05:7.0
data(:,1) = 0.0001:0.0001:0.9999; % just following Brent's files
data(:,2) = ro.*nthroot(((1-data(:,1)).^2)./data(:,1),6);  % Estimated FRET from Vogel Methods 2012 and Vogel PLoS ONE 2012
filename = sprintf('e_from_rc_ro_%0.6f.mat', ro);
save(filename, 'data'); 
end