function output = FRET_vs_Na(na_min,na_step,na_max,ro_min,ro_step,ro_max,rfp_min,rfp_step,rfp_max,lp_min,lp_step,lp_max,n)
na_vector = na_min:na_step:na_max;
ro_vector = ro_min:ro_step:ro_max;
rfp_vector = rfp_min:rfp_step:rfp_max;
lp_vector = lp_min:lp_step:lp_max;

%% FUNCTION DESCRIPTION
% This wrapper function calls and compiles the output from
% construct_design.m over all input linker lengths, forster
% distances, FP radii, and persistence lengths, n-times.
% Ultimately, the output of this function describes the mechanical response
% of molecular tension sensors in unloaded conditions
%
%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%
% na_vector  -> make sure this is 1) a vector 2) in units of amino acids
% ro_vector  -> make sure this is 1) a vector 2) in units of nm
% rfp_vector -> make sure this is 1) a vector 2) in units of nm
% lp_vector  -> make sure this is 1) a vector 2) in units of nm
% n          -> 1e3 is an absolute minimum (quick and dirty); 1e5 is very smooth (more is probably overkill)(usually do 1e4)
%
%%%%%%%%%%%%%%%%%%%%% OUTPUTS %%%%%%%%%%%%%%%%%%%%%
% SEE construct_design.m DESCRIPTION FOR FURTHER DETAILS
% Saves input parameters that correpond to individual unloaded
% linker mechanical responses:
%     output.ro:       Förster distance(nm)
%     output.rf:       Fluorescent protein radius(nm)
%     output.lp:       Persistence length(nm)
%     output.na:       Number of amino acids in polypeptide linker(amino acids)
% Saves variables that describe the linker mechanical response at F = 0 pN
%     output.re       Polypeptide end-to-end distance (nm)
%     output.rc1-4    Chromophore separation distances (nm)
%                     Multiple outputs because this calc depends on:
%                     (1) re_2_rc conversion assumption
%                     (2) Numerical averaging
%     output.e1-4D/S  FRET efficiencies (decimal) - rows in order of
%                     increasing linker length; columns in order
%                     corresponding to parameter vectors.
%                     Multiple outputs because this calc depends on:
%                     (1) re_2_rc conversion assumption
%                     (2) Numerical averaging
%                     (3) Assumption as to whether you are in the "dynamic
%                     isotropic regime" of FRET calculation or the "static
%                     isotropic regime" of FRET

%% Add lookup table folders to path
if exist('FRET-NA Model Results','dir')~=7
    mkdir('FRET-NA Model Results')
end
addpath(genpath(pwd));

%% Define constant parameters
cn = 1;
temp = 37;
fsteps = 1;
fjc = 0;
wlc = 1;
cyl = 0;
bo = 0.38;
conlc = 0;
fmax = 0;
fmin = 0;
zfrc = 0;
k = 0;
z = 0;
iterations = length(ro_vector).*length(rfp_vector).*length(lp_vector);
iterations1 = iterations.*length(na_vector);
h = waitbar(0,'Please wait...');

%% Pre-allocate space for individual outputs
output.ro = zeros(1, iterations);
output.rfp = zeros(1, iterations);
output.lp = zeros(1, iterations);
output.na = zeros(1, iterations);
output.re = zeros(length(na_vector), iterations);
output.rc1 = zeros(length(na_vector), iterations);
output.rc2 = zeros(length(na_vector), iterations);
output.rc3 = zeros(length(na_vector), iterations);
output.rc4 = zeros(length(na_vector), iterations);
output.e1D = zeros(length(na_vector), iterations);
output.e1S = zeros(length(na_vector), iterations);
output.e2D = zeros(length(na_vector), iterations);
output.e2S = zeros(length(na_vector), iterations);
output.e3D = zeros(length(na_vector), iterations);
output.e3S = zeros(length(na_vector), iterations);
output.e4D = zeros(length(na_vector), iterations);
output.e4S = zeros(length(na_vector), iterations);

%% Loop over ro, rfp, lp, na (single force)
for a = 1:length(ro_vector)
    ro = ro_vector(a);
    rofn = sprintf('e_from_rc_ro_%0.6f.mat', ro);
    for b = 1:length(rfp_vector)
        rfp = rfp_vector(b);
        rfpfn = sprintf('quad_lc_%0.6f_%0.6f.mat', rfp, rfp);
        for c = 1:length(lp_vector)
            lp = lp_vector(c);
            k = k+1;
            output.ro(k) = ro;
            output.rfp(k) = rfp;
            output.lp(k) = lp;
            
            for d = 1:length(na_vector)
                na = na_vector(d);
                output.na(d,k) = na;
                z = z+1;
                
                % Calculations always only at F = 0
                res = construct_design(cn, na, lp, rfp, rfpfn, rofn, ro, n, temp, fsteps,...
                    fjc, wlc, cyl, bo, conlc, fmax, fmin, zfrc);
                
                % Assemble individual outputs
                output.re(d,k) = res(:,3);
                output.rc1(d,k) = res(:,4);
                output.rc2(d,k) = res(:,5);
                output.rc3(d,k) = res(:,6);
                output.rc4(d,k) = res(:,7);
                output.e1D(d,k) = res(:,8);
                output.e1S(d,k) = res(:,9);
                output.e2D(d,k) = res(:,10);
                output.e2S(d,k) = res(:,11);
                output.e3D(d,k) = res(:,12);
                output.e3S(d,k) = res(:,13);
                output.e4D(d,k) = res(:,14);
                output.e4S(d,k) = res(:,15);
                
                waitbar(z./iterations1,h);
            end
        end
    end
end
%% Construct Complicated Savename
na_min_str = strrep(num2str(na_min),'.','p');
na_step_str = strrep(num2str(na_step),'.','p');
na_max_str = strrep(num2str(na_max),'.','p');
if strcmpi(na_min_str,na_max_str)
    na_str = strcat('NA',na_min_str,'AAs_');
else
    na_str = strcat('NA',na_min_str,'by',na_step_str,'to',na_max_str,'AAs_');
end

ro_min_str = strrep(num2str(ro_min),'.','p');
ro_step_str = strrep(num2str(ro_step),'.','p');
ro_max_str = strrep(num2str(ro_max),'.','p');
if strcmpi(ro_min_str,ro_max_str)
    ro_str = strcat('Ro',ro_min_str,'nm_');
else
    ro_str = strcat('Ro',ro_min_str,'by',ro_step_str,'to',ro_max_str,'nm_');
end

rfp_min_str = strrep(num2str(rfp_min),'.','p');
rfp_step_str = strrep(num2str(rfp_step),'.','p');
rfp_max_str = strrep(num2str(rfp_max),'.','p');
if strcmpi(rfp_min_str,rfp_max_str)
    rfp_str = strcat('Rfp',rfp_min_str,'nm_');
else
    rfp_str = strcat('Rfp',rfp_min_str,'by',rfp_step_str,'to',rfp_max_str,'nm_');
end

lp_min_str = strrep(num2str(lp_min),'.','p');
lp_step_str = strrep(num2str(lp_step),'.','p');
lp_max_str = strrep(num2str(lp_max),'.','p');
if strcmpi(lp_min_str,lp_max_str)
    lp_str = strcat('Lp',lp_min_str,'nm_');
else
    lp_str = strcat('Lp',lp_min_str,'by',lp_step_str,'to',lp_max_str,'nm_');
end

n_str = strrep(num2str(n),'.','p');
savename = strcat('FRET_vs_Na_',na_str,ro_str,rfp_str,lp_str, n_str, 'iterations');
save(fullfile(pwd,'FRET-NA Model Results',savename),'output');
close
end