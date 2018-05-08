function res = construct_design(cn, na, lp, rfp, rfpfn, rofn, ro, n, temp, fsteps,...
    fjc, wlc, cyl, bo, conlc, fmax, fmin, zfrc)

%% FUNCTION DESCRIPTION
%  This function predicts both the means and distributions of the following
%  as a function of applied force:
%       fz - Force "felt" by the linker (input)
%       rz - Polypeptide linker extension (output)
%       re - Polypeptide linker end-to-end distance (output)
%       rc - Fluorescent protein chromophore separation distance (multiple outputs, 4)
%       e  - FRET efficiency (multiple outputs, 8)
%       f  - Force applied (multiple outputs since critical force scale
%            depends on the four different rcs)

%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%
% Note: if variable = 1, model type is in use, otherwise, it is not
% Note: all lengths in nm
%       cn        - Parameter utilized in the context of the Gaussian and
%                   FJC models to represent the mechanical behavior of a
%                   real chain, or persistence length for WLC
%       na        - Number of amino acids (or Kuhn lengths in the FJC)
%       lp        - Persistence lengt
%       rfp       - Fluorophore radius
%       rfpfn     - Name of the file containing conversion of end to end distance of
%                   polypeptide (re) to chromophore separation distance (rc)
%       rofn      - Name of the file containing conversion of chromophore
%                   separation distance (rc) to FRET efficiency (e) in "static
%                   isotropic regime"
%       ro        - Forster distance
%       n         - number of points in distributions
%       temp      - temperature in degrees C
%       fsteps    - number of force steps
%       fjc = 1   - use Freely-Jointed Chain model instead of Gaussian Chain model (default)
%       wlc = 1   - use Worm-Like Chain model instead of Gaussian Chain model (default)
%       cyl = 1   - use Cylindrical model instead of Gaussian Chain model (default)
%       bo        - Monomer size for gaussian or WLC, Kuhn for FJC
%                   aka, length of amino acid in nm (0.71 nm for DNA)
%       conlc = 1 - constrain to contour length
%       fmax      - maximal force
%       fmin      - minimum force
%       zfrc      - chromophore separation distance (rc) at F = 0;

%%%%%%%%%%%%%%%%%%%%% OUTPUTS %%%%%%%%%%%%%%%%%%%%%
% Only one potential average polymer end to end distance (re) and
% normalized extension (rz)
%    VARIABLE     Description
%       re     -  <p(re)>
%       rz     -  integral[p(re)] |-lc to lc
%
% One can get a variety of chromophore separation values depending on
% (A) Which re to rc conversion technique is utilized
%        Brenner: rc = re + C
%        Hoffman: rc = heuristic sum of re + rfp
% (B) Whether numerical averages <re> or distributions p(re) are utilized
%    VARIABLE    Avg./Distrib.    Re2Rc Conv.
%      rc1     -    avg            re + C
%      rc2     -    avg         heuristic sum
%      rc3     -    dist           re + C
%      rc4     -    dist        heuristic sum
%
% From these rc 1 thru 4 values, 8 FRET efficiencies can be calculated,
% depending on whether you assume you are in the dynamic isotropic regime
% (classical assuption, use classic Forster equation) or you are in the
% static isotropic regime (more recent Monte Carlo simultations indicate
% this might be more appropriate for bulky FPs, Vogel, et al., Methods 2012)
%    VARIABLE    Avg./Distrib.    Re2Rc Conv.   FRET regime
%      e1D     -    avg            re + C         dynamic
%      e1S     -    avg            re + C         static
%      e2D     -    avg         heuristic sum     dynamic
%      e2S     -    avg         heuristic sum     static
%      e3D     -    dist           re + C         dynamic
%      e3S     -    dist           re + C         static
%      e4D     -    dist        heuristic sum     dynamic
%      e4S     -    dist        heuristic sum     static
%
% Again, from these rc 1 thru 4 values, 4 corresponding Force_applieds can
% be calculated
%       f1     -   Fapplied based on Fcrit from rc1
%       f2     -   Fapplied based on Fcrit from rc2
%       f3     -   Fapplied based on Fcrit from rc3
%       f4     -   Fapplied based on Fcrit from rc4

%% Add pertinent folders to path and pre-allocate space for outputs
addpath('Support Functions');
f = zeros(1,fsteps);
rz = zeros(1,fsteps);
re = zeros(1,fsteps);
rc1 = zeros(1,fsteps);
rc2 = zeros(1,fsteps);
rc3 = zeros(1,fsteps);
rc4 = zeros(1,fsteps);
e1D = zeros(1,fsteps);
e1S = zeros(1,fsteps);
e2D = zeros(1,fsteps);
e2S = zeros(1,fsteps);
e3D = zeros(1,fsteps);
e3S = zeros(1,fsteps);
e4D = zeros(1,fsteps);
e4S = zeros(1,fsteps);

%% Define polymer model
if fjc == 1
    type = 1;
elseif wlc == 1
    type = 2;
elseif cyl == 1
    type = 3;
else
    type = 0;
    disp('Gaussian Chain model selected by default');
end

%% Define Additional Parameters
% Force vector (this is the force "felt" by the linker, Fz)
fz = linspace(fmin, fmax, fsteps);

% Thermal energy and critical force scales
kbt = (temp+273.15).*1.3806503.*10.^(-23).*10.^21;
if zfrc == 0
    F_crit = 0;
else
    F_crit = zeros(1,4);
    for j = 1:4
        F_crit(j) = kbt./zfrc(j);
    end
end

% Contour lengths
lc = bo.*na;
switch type
    case 0
        lmax = na.*bo.*4;
    case 1
        lmax = 0.99.*lc;
    case 2
        lmax = 0.99.*lc;
    case 3
        lmax = 0.99.*lc;
end

%% Load conversion lookup tables
% Used in re_2_rc conversion (based on input rfp)
table1 = load(rfpfn);
conv_rfp = table1.data;

% Used in FRETE_static conversion (based on input ro)
table2 = load(rofn);
conv_ro = table2.data;

%% Get the distributions at zero force
fzero = 0;
% Getting P(re)
noend = 0; e2e = 1;
r = linspace(0, lmax, 1e3);
switch type
    case 0
        rfunc = gc_func(r, na, cn, bo, fzero, temp, noend);
    case 1
        rfunc = fjc_func(r, na, bo, fzero, temp, noend);
    case 2
        rfunc = wlc_becker_func(r, lp, lc, e2e);
    case 3
        rfunc = cylinder_func(r, lp, lc, e2e);
end

% Getting P(R)
noend = 1; e2e = 0;
rneg = linspace(-lmax, lmax, 2e3);
switch type
    case 0
        efunc = gc_func(rneg, na, cn, bo, fzero, temp, noend);
    case 1
        efunc = fjc_func(rneg, na, bo, fzero, temp, noend);
    case 2
        efunc = wlc_becker_func(rneg, lp, lc, e2e);
    case 3
        efunc = cylinder_func(rneg, lp, lc, e2e);
end

%% Loop over force
for i = 1:fsteps
    % Calculate the extension, rz function,  -lc -> lc
    fefunc = efunc.*exp(fz(i).*rneg./kbt);
    rz(i) = ext_ch(rneg, fefunc, lc, type, conlc); % older method, ok to use because we never use P(rz) in downstream calculations
    
    % Calcualte end-to-end distance, re function,  0 -> lc
    frfunc = rfunc.*exp(fz(i).*r./kbt);
    normf = trapz(r, frfunc);
    frfunc = frfunc./normf;
    pre = ch_dist(n, r, frfunc, lc, type, conlc); % Take distribution(re) and make it a list(re) using inverse transform sampling method
    in = pre ~= Inf; % Remove infinities
    pre = pre(in); % p(re)
    re(i) = mean(pre); % re
    
    % Convert list(re) to list(rc)--chromophore to chromophore distance
    % VARIABLE   CALCULATION               Avg./Distrib.    Re2Rc Conv.
    rc1(i)   =  re(i) + 2*rfp;             %    avg             re + C
    rc2(i)   =  re_2_rc(re(i),conv_rfp);   %    avg          heuristic sum
    prc3  =  pre + 2*rfp;                  %    dist            re + C
    prc4  =  re_2_rc(pre, conv_rfp);       %    dist         heuristic sum
    % Avgs of distributions
    rc3(i) = mean(prc3);
    rc4(i) = mean(prc4);
    
    % Calculate FRET efficiencies
    % VARIABLE     CALCULATION           Avg./Distrib.     Re2Rc Conv.    FRET regime
    e1D(i)  = FRETE_dynamic(rc1(i), ro);     %avg            re + C         dynamic
    e1S(i)  = FRETE_static(rc1(i),conv_ro);  %avg            re + C         static
    e2D(i)  = FRETE_dynamic(rc2(i), ro);     %avg         heuristic sum     dynamic
    e2S(i)  = FRETE_static(rc2(i),conv_ro);  %avg         heuristic sum     static
    pe3D = FRETE_dynamic(prc3, ro);          %dist           re + C         dynamic
    pe3S = FRETE_static(prc3,conv_ro);       %dist           re + C         static
    pe4D = FRETE_dynamic(prc4, ro);          %dist        heuristic sum     dynamic
    pe4S = FRETE_static(prc4,conv_ro);       %dist        heuristic sum     static
    % Avgs of distributions
    e3D(i) = mean(pe3D);
    e3S(i) = mean(pe3S);
    e4D(i) = mean(pe4D);
    e4S(i) = mean(pe4S);
    
    % Calculate force applied (F) from the force "felt" by the linker (Fz)
    if length(F_crit) > 1
        if fz(i) < 0
            for j = 1:4
                fz(i) = -fz(i); %% makes it positive
                f(j,i) = ((fz(i) + F_crit(j)).^2 -F_crit(j).^2).^(1./2); %% calculates as if it were positive
                fz(i) = -fz(i); %% changes the felt force back
                f(j,i) = -f(j,i); %% changes the applied force to it's actual value
            end
        else
            for j = 1:4
                f(j,i) = ((fz(i) + F_crit(j)).^2 -F_crit(j).^2).^(1./2);
            end
        end
    else
        if fz(i) < 0
            fz(i) = -fz(i); %% makes it positive
            f(i) = ((fz(i) + F_crit).^2 -F_crit.^2).^(1./2); %% calculates as if it were positive
            fz(i) = -fz(i); %% changes the felt force back
            f(i) = -f(i); %% changes the applied force to it's actual value
        else
            f(i) = ((fz(i) + F_crit).^2 -F_crit.^2).^(1./2);
        end
    end
end

%% Compile output
if length(F_crit)>1
    res = [fz', rz', re', rc1', rc2', rc3', rc4', e1D', e1S', e2D', e2S', e3D', e3S', e4D', e4S', f(1,:)', f(2,:)', f(3,:)', f(4,:)'];
else
    res = [fz', rz', re', rc1', rc2', rc3', rc4', e1D', e1S', e2D', e2S', e3D', e3S', e4D', e4S', f'];
end
end