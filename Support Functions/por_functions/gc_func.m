%% This function describes the mechanics of polymers using the Gaussian Chain Model
% Uses the equation for Gaussian Chain extension at zero force in z direction
% inputs: r is a vector of end-to-end distances
%         lp is the persistence length
%         lc is the contour length
% outputs:
%         if e2e is 0, then gor = P(R) = end-to-end VECTOR probability density function
%         if e2e is 1, then gor = P(re) = end-to-end DISTANCE probability density function

%% Function gc_func
function pxo = gc_func(r, na, cn, bo, f, temp, noend)
kbt = (temp + 273.15).*1.3806503.*10.^(-23).*10.^21;

% Determine either P(R) or P(re)
if noend == 1 % P(R) = end-to-end VECTOR probability density function
    pxo = (3./(2.*pi.*cn.*na.*bo.^2)).^(3./2).*exp(-3.*r.^2./(2.*cn.*na.*bo.^2)); 
%     pxo = (3./(2.*pi.*cn.*na.*bo.^2)).^(3./2).*exp(-3./(2.*cn.*na.*bo.^2).*(r-(cn.*na.*bo.^2.*f)./(3.*kbt)).^2);
else % P(re) = end-to-end DISTANCE probability density function
    pxo = (4.*pi.*r.^2).*(3./(2.*pi.*cn.*na.*bo.^2)).^(3./2).*exp(-3.*r.^2./(2.*cn.*na.*bo.^2)); 
%     pxo = (4.*pi.*r.^2).*(3./(2.*pi.*cn.*na.*bo.^2)).^(3./2).*exp(-3./(2.*cn.*na.*bo.^2).*(r-(cn.*na.*bo.^2.*f)./(3.*kbt)).^2);
end

% Normalize area under curve to 1
normf = trapz(r, pxo);
pxo = pxo./normf;