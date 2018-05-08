%% Function for the second Daniels approximation of the WLC
% Option 2 for soft polymers and/or small extensions lc/lp > 8
% DANIELS DOESN'T WORK FOR FORCES NEAR FULL EXTENSION
% This is how it is in Brent's code--verified from other papers (Zhou, Yamakawa)
% inputs: r is a vector of end-to-end distances
%         lp is the persistence length
%         lc is the contour length
% outputs:
%         if e2e is 0, then gor = P(R) = end-to-end VECTOR probability density function
%         if e2e is 1, then gor = P(re) = end-to-end DISTANCE probability density function
%% Function wlc_daniels_func
function gor = wlc_daniels_func(r, lp, lc, f, temp, noend, noise)
% Constants
rn = r./lc; % r is scaled to lc
nele = numel(r);
gor = zeros(1,nele);
kbt = (temp+273.15).*1.3806503.*10.^(-23).*10.^21;

% Eight terms in second Daniels approx
%   Different authors did this differently, I'm 90% sure that wilhelm did it
%   correctly and that all the other authors may have used the correct
%   equations, but published incorrect versions of them
author = 'Wilhelm';
if strcmpi(author,'Zhou') || strcmpi(author,'Yamakawa') || strcmpi(author,'Hoffman') % questionable signs
    t1 = (-5.*lp)./(4.*lc);
    t2 = (2.*r.^2)./(lc.^2);
    t3 = (-33.*r.^4)./(80.*lp.*lc.^3);
    t4 = (-79.*lp.^2)./(160.*lc.^2);
    t5 = (-329.*r.^2.*lp)./(120.*lc.^3);
    t6 = (6799.*r.^4)./(1600.*lc.^4);
    t7 = (-3441.*r.^6)./(2800.*lc.^5.*lp);
    t8 = (1089.*r.^8)./(12800.*lp.^2.*lc.^6);
elseif strcmpi(author,'LaCroix') || strcmpi(author,'Wilhelm') % correct signs
    t1 = (5.*lp)./(4.*lc);
    t2 = (2.*r.^2)./(lc.^2);
    t3 = (-33.*r.^4)./(80.*lp.*lc.^3);
    t4 = (-79.*lp.^2)./(160.*lc.^2);
    t5 = (-329.*r.^2.*lp)./(120.*lc.^3);
    t6 = (6799.*r.^4)./(1600.*lc.^4);
    t7 = (-3441.*r.^6)./(2800.*lc.^5.*lp);
    t8 = (1089.*r.^8)./(12800.*lp.^2.*lc.^6);
end

% Compile terms
corterm = 1+t1+t2+t3+t4+t5+t6+t7+t8;
gor = ((3)./(4.*pi.*lc.*lp)).^(1.5).*exp((-3.*r.^2)./(4.*lp.*lc)).*corterm;

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