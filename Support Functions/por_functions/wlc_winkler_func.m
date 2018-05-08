%% This function describes the mechanics of polymers using some mystery approximation of the WLC Model that Brent gave to Drew
% Comes from Winkler et al., J Chem Phys 2003
% Conclusion as of 8/7/17 was that this may have been useful before Becker
% et al., 2010 paper came out - this use to be called "3rd magic function
% from Brent"
% Option 1 for soft polymers and/or small extensions lc/lp > 8
% inputs: r is a vector of end-to-end distances
%         lp is the persistence length
%         lc is the contour length
% outputs:
%         if e2e is 0, then gor = P(R) = end-to-end VECTOR probability density function
%         if e2e is 1, then gor = P(re) = end-to-end DISTANCE probability density function

%% Function wlc_magic_func
function gor = wlc_winkler_func(r, lp, lc, f, temp, noend, noise)
% Constants
flag = 1;
rn = r./lc; % r is scaled to lc
nele = numel(r);
nterms = 10;
gor = zeros(1,nele);
kbt = 1;
%kbt = (temp+273.15).*1.3806503.*10.^(-23).*10.^21;

% 3rd magic function from Brent
gor = ((1-rn.^2).^(-1.5)).*((2-rn.^2).^(-3)).*exp(-3.*lc./4./lp./(1-rn.^2));
for i = 1:numel(gor)
    if gor(i)./Inf <= 0
        gor = gor;
    else
        gor(i) = 0;
    end
end

% Determine either P(R) or P(re)
if noend == 1
    gor = gor; % P(R) = end-to-end VECTOR probability density function
else
    gor = 4.*pi.*rn.^2.*gor; % P(re) = end-to-end DISTANCE probability density function
end

% Normalize area under curve to 1
normf = trapz(rn, gor);
gor = gor./normf;
end