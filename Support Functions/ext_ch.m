%% older method for force extension curve--based on taking integral of P(r,f)
function x = ext_ch(r, efunc, lc, type, conlc)

% Never have to do this...
if conlc == 1
    wglc = find(r >= lc);
    if numel(wglc) >= 1
        tot = sum(efunc(wglc));
        efunc(wglc(1)) = efunc(wglc(1)) + tot;
        efunc(wglc(2:end)) = 0;
    end
end

% Normalize [P(re) from -Lc to Lc] to 1
normf = trapz(r, efunc);
efunc = efunc./normf;

%
x = trapz(r, r.*efunc);
end