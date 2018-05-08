%% This function describes the mechanics of Worm-like Chains based on Becker, Rosa, and Everaers (2010)
% inputs: r is a vector of end-to-end distances
%         lp is the persistence length
%         lc is the contour length
% outputs: noend = 0; e2e = 1;
%         if e2e is 0, then gor = P(R) = end-to-end VECTOR probability density function
%         if e2e is 1, then gor = P(re) = end-to-end DISTANCE probability density function

%% Function wlc_becker_func
function gor = wlc_becker_func(r, lp, lc, e2e)  % gor is detailed above
r    = r./lc;                                   % normalized length (0 to 1)
k    = lp./lc;                                  % normalized stiffness term
c    = c_func(lp, lc);                          % c and d are free parameters dependent on the stiffness
Qa   = Qa_func(r, lp, lc);                      % Qa is responsible for retaining the correct
                                                % assymtotes at r->1 and r->0
Qb   = Qb_func(r, lp, lc);                      % Qb accounts for the exponential suppresion of high
                                                % energy closed-loop configurations of rigid rods
lead = (1-c.*r.^2).^(5/2);                      % the leading term tunes the cross-over region
                                                % at intermediate stiffnesses
Q    = lead.*Qa.*Qb;

% Determine either P(R) or P(re)
if e2e == 0
    gor = Q; % P(R) = end-to-end VECTOR probability density function
else
    gor = 4.*pi.*r.^2.*Q; % P(re) = end-to-end DISTANCE probability density function
end

% Make sure no NaNs exist
for z = 1:length(r)                          
    if (gor(z)./Inf) == 0
        gor(z) = gor(z);
    else
        gor(z) = 0;
    end
end

% Normalize area under curve to 1
normf = trapz(r, gor);
gor = gor./normf;
end

%% Subfunction Qa--responsible for retaining the correct assymtotes at the boundaries
function Qa = Qa_func(r, lp, lc)
k   = lp/lc;
% cij = [-3/4 23/64 -7/64 -1/2 7/16 -9/16]; % BAD terms derived from matching orders of r and k
cij = [-3/4 23/64 -7/64 -1/2 17/16 -9/16]; % FIXED terms derived from matching orders of r and k

term1 = Jsyd_func(lp, lc);                 % the Stockmayer J-factor
term2 = (1-r.^2).^(-5/2);                  % a tuning factor

n     = 0;                                 % estabilishing a counter for the loop
term3 = 0;                                 % establishing the variable to be modified

for i = -1:0
    for m = 1:3
        n = n + 1;
        term3 = term3 + (cij(n).*k.^(i).*r.^(2.*m))./(1-r.^2);
    end
end
term3 = exp(term3);                        % the main term for Qa

Qa    = term1.*term2.*term3;               % combining all terms
end

%% Subfunction Qb--accounting for the exponential suppression of high energy closed-loop configurations of rigid rods
function Qb = Qb_func(r, lp, lc)
k = lp/lc;
a = 14.054;                                               % related to the classical closed-loop energy
b = 0.473;                                                % related to the classical closed-loop energy
d = d_func(lp, lc);                                       % free parameter dependent on stiffness k

Io = -(d.*k.*a.*(1 + b).*r)./(1 - b.^2.*r.^2);            % this term feeds into the modified
                                                          % bessel function of the first kind

bessel = besseli(0, Io);                                  % evaluates the bessel function
term   = exp(-(d.*k.*a.*b.*(1+b).*r.^2)./(1-b.^2.*r.^2)); % the main term for Qb

Qb     = term.*bessel;                                    % combining all terms
end

%% Subfunction Jsyd-- calculates teh Stockmayer J-factor used in Qa calculation
function Jsyd = Jsyd_func(lp, lc) % function for the Stockmayer J-factor calculated by Shimada and Yamakawa
k = lp/lc;
a = 14.054;                       % related to the classical closed-loop energy

if k > (1/8)                      % non-soft polymers
%     Jsyd = 112.04.*k.^2.*exp(0.246/(k-a.*k)); % BAD
    Jsyd = 112.04.*k.^2.*exp((0.246/k)-a.*k); % FIXED
else                              % soft polymers
    Jsyd = (3/(4.*pi.*k)).^(3/2).*(1-5.*k/4);
end

end

%% Subfunction c--calculates the free parameter c
function c = c_func(lp, lc) % function for the free parameter c
k = lp/lc; % lp/lc is a normalized effective stiffness

c = 1 - (1 + (0.38.*k.^(-0.95)).^(-5)).^(-1/5);
end

%% Subfunction d--calculates the free parameter d
function d = d_func(lp, lc) % function for the free parameter d
k = lp/lc;

if k < (1/8)                % exception for soft polymers
    d = 1;
else                        % all other cases
    d = 1 - ((0.177/(k-0.111)) + 6.40.*(k-0.111).^(0.783))^(-1);
end
end