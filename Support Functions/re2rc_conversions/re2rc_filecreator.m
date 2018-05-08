%% Script to make re -> rc conversion files
clear;
for r1 = 0.05:0.05:4.95
data(:,1) = linspace(0, 100*49.98, 100*2500); % just following Brent's files
data(:,2) = sqrt((r1+r1).^2 + data(:,1).^2);  % quadrature formula
filename = sprintf('quad_lc_%f_%f.mat', r1, r1);
save(filename, 'data'); 
end