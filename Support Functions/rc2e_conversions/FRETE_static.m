%% Function to Convert End to End Distance to FRET Efficiency [aka P(rc) to P(E)]
% Uses the lookuptable that describes the relationship between (r/ro)^6 and
% FRET efficiency in the static random isotropic averaging regime
function E = FRETE_static(x, ro_conv)
x(x>max(ro_conv(:,2))) = max(ro_conv(:,2));
x(x<min(ro_conv(:,2))) = min(ro_conv(:,2));
E = interp1(ro_conv(:,2), ro_conv(:,1), x);
end