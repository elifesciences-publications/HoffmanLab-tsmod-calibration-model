%% This function describes the mechanics of polymers using the Wilhem-Frey series expansion of the WLC
% works well at high extensions (not at lower ones) and high stiffness, lc/lp <= 8
% inputs: r is a vector of end-to-end distances
%         lp is the persistence length
%         lc is the contour length
% outputs:
%         if e2e is 0, then gor = P(R) = end-to-end VECTOR probability density function
%         if e2e is 1, then gor = P(re) = end-to-end DISTANCE probability density function

%% Function wlc_WF_func
function gor = wlc_WF_func(r, lp, lc, f, temp, noend, noise)
% Constants
flag = 1;
rn = r./lc; % r is scaled to lc
nele = numel(r);
nterms = 10;
gor = zeros(1,nele);
kbt = 1;
%kbt = (temp+273.15).*1.3806503.*10.^(-23).*10.^21;

% Wilhelm Result
m = 1:nterms;
for i = 1:(nele-1)
    thres = (1-.2.*lc./lp);
    if (rn(i) < thres) && (flag == 1)
        for m = 1:nterms
            %how it is in brent's code:
            term = (pi.^2).*(m.^2).*(-1).^(m+1).*exp(-kbt.*lp./lc.*(pi.^2).*(m.^2).*(1-rn(i)));
            term = 2.*kbt.*lp./lc./4./pi.*term;
            gor(i) = gor(i) + term;
        end
    else
        if (flag == 1)
            pgor = pgor(lp, lc, rn, m, i-1, nterms); % gets previous g(t)
            if i == 1
                tnormf = 1;
                flag = 0;
            else
                tnormf = gor(i-1)./pgor;
                flag = 0;
            end
        end
        for m = 1:nterms
            term = lr_approx(lp, lc, rn, m, i, nterms);
            term = term.*tnormf;
            gor(i) = gor(i) + term;
        end
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

%% Subfunction lr_approx
function term = lr_approx(lp, lc, rn, m, i, nterms)
if i == 0
    grp = (m-0.5)./sqrt(lp./lc);
    herm = hermite(grp);
    term = (1./(lp./lc).^(1.5)).*exp(-(m-0.5).^(2).*lc./lp).*herm;
    term = lp./lc./8./pi./sqrt(pi).*term;
else
    grp = (m-0.5)./sqrt(lp./lc.*(1-rn(i)));
    herm = hermite(grp);
    term = (1./(lp./lc.*(1-rn(i))).^(1.5)).*exp(-(m-0.5).^(2).*lc./lp./(1-rn(i))).*herm;
    term = lp./lc./8./pi./sqrt(pi).*term;
end
end

%% Subfunction hermite
function herm = hermite(x)
herm = 4.*(x.^2) -2;
end

%% Subfunction pgor (previous gor)
function gor = pgor(lp, lc, rn, m, i, nterms)
gor = 0;
for m = 1:nterms
    term = lr_approx(lp,lc,rn,m,i, nterms);
    gor = gor + term;
end
end