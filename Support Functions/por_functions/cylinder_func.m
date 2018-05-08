%% This function describes the mechanics of Worm-like Chains based on some kind of simple cylindrical or analagous mechanics
% Tested only for compression, but didn't see much difference between this
% and the negative version of the Becker et al., 2010 model.
% inputs: r is a vector of end-to-end distances
%         lp is the persistence length
%         lc is the contour length
% outputs:
%         if e2e is 0, then gor = P(R) = end-to-end VECTOR probability density function
%         if e2e is 1, then gor = P(re) = end-to-end DISTANCE probability density function

%% Function cylinder_func
function gor = cylinder_func(r, lp, lc, e2e)

r = r./lc;
y = lp./lc./2.*(1-r.^2);

% Determine either P(R) or P(re)
if e2e == 0
    gor = exp(-pi.^2.*y - 1./pi./y); % P(R) = end-to-end VECTOR probability density function
else
    gor = 4.*pi.*r.^2.*exp(-pi.^2.*y - 1./pi./y); % P(re) = end-to-end DISTANCE probability density function
end

% Normalize area under curve to 1
normf = trapz(r, gor);
gor   = gor./normf;
end