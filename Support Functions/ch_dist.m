%% ch_dist function
% turns a distribution into a randomly selected set of numbers
% weighted by the distribution (inverse transform sampling method?)
function y = ch_dist(n,r,cfunc,lc,type,conlc)
if type==5
%     cfunc   = r.*cfunc;
    cdf     = cumsum(cfunc);              % r(2) normalizes step size
    cdf     = cdf./max(cdf);
    end_cdf = find(cdf==1);
    cdf     = cdf(1:end_cdf(1));          % cut off the cdf once it hits 1
    start_cdf = find(cdf==0);
    cdf = cdf(start_cdf(end):end);        % cut off the cdf once it hits 0 too
    r       = r(1:end_cdf(1));            % r needs to be the same size as the cdf
    r       = r(start_cdf(end):end);      % r needs to be the same size as the cdf
    u       = rand(1,n);                  % random uniform distribution of 10,000 numbers
    y       = zeros(1, n);
else
    cdf     = cumsum(cfunc);       % r(2) normalizes step size
    cdf     = cdf./max(cdf);
    end_cdf = find(cdf==1);
    cdf     = cdf(1:end_cdf(1));   % cut off the cdf once it hits 1
    r       = r(1:end_cdf(1));     % r needs to be the same size as the cdf
    u       = rand(1,n);           % random uniform distribution of 10,000 numbers
    y       = zeros(1, n);
end

% Shouldn't ever have to do this due to above calculations by definition
if max(cdf) > 1
    w = find(cdf <= 1);
    cdf = cdf(w);
    r = r(w);
    y = y(w);
end

% 10,000 times, calculate Y(
for i = 1:n
    y(i) = interp1(cdf, r, u(i));
    if conlc == 1
        y(y>lc) = lc;
    end
end
y(y>lc) = lc;

%% remove NaN
real = y./Inf==0;
y = y(real);

end
