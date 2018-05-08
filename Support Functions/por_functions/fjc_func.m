%% This function describes the mechanics of polymers using the Freely-Jointed Chain Model
% Uses the equation for approximation of FJC extension at zero force in z direction
% inputs: r is a vector of end-to-end distances
%         lp is the persistence length
%         lc is the contour length
% outputs:
%         if e2e is 0, then gor = P(R) = end-to-end VECTOR probability density function
%         if e2e is 1, then gor = P(re) = end-to-end DISTANCE probability density function

%% Function fjc_func
function pxo = fjc_func(r, na, bo, f, temp, noend)

kbt = (temp + 273.15).*1.3806503.*10.^(-23).*10.^21;
beta = inv_lang(r, na, bo);

loggoo = -na.*(r./na./(bo).*beta+log(beta./sinh(beta)));
goo    = exp(loggoo);

% Determine either P(R) or P(re)
if noend == 1
    pxo = goo; % P(R) = end-to-end VECTOR probability density function
    pxo(1) = 0;
    pxo(end) = 0;
    %exp(f.*r./kbt).*goo;
else
    pxo = (4.*pi.*r.^2).*goo; % P(re) = end-to-end DISTANCE probability density function
    %(4.*pi.*r.^2).*exp(f.*r./kbt).*goo;
end

% finding the non-infinite max
CHECK = zeros(1, numel(r));
for a = 1:length(r)
    if pxo(a) < Inf
        CHECK(a) = pxo(a);
    end
end
real_max = max(CHECK);

% replacing infinite/NaN numbers
for b = 1:length(r)
    if pxo(b) == Inf
        pxo(b) = 0;
    elseif pxo(b) > -1
        pxo(b) = pxo(b);
    else
        pxo(b) = real_max;
    end
    
    % if pxo(b) < Inf
    %     pxo(b) = pxo(b);
    % else
    %     pxo(b) = 0;
    % end
end

% Normalize area under curve to 1
normf = trapz(r, pxo);
pxo = pxo./normf;
end

%% Subfunction inv_lang
function res = inv_lang(r, na, bo)
nele = 3e3;
nele_vect = linspace(0,nele-1,nele);
x = (1.01.^nele_vect)./1e2; % set of log-spaced numbers
xneg = -x;
x = [xneg x];
lang = 1./tanh(x) - 1./x;
vec = r./(na.*bo);
disp(max(vec))
if max(vec) > 1
    warning('renormalize the distance')
end
res = r;
res(:) = 0;
for i = 1:length(r)
    if vec(i) > max(lang)
        res(i) = max(lang);
    else
        res(i) = interp1(lang, x, vec(i));
    end
end
end