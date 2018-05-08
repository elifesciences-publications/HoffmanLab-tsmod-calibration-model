%% Function to Convert FRET Efficiency to End to End Distance [aka P(E) to P(rc)]
% Uses the lookuptable that describes the relationship between (r/ro)^6 and
% FRET efficiency in the static random isotropic averaging regime
function x = inv_FRETE_static(E, ro_conv)
E(E>max(ro_conv(:,1))) = max(ro_conv(:,1));
E(E<min(ro_conv(:,1))) = min(ro_conv(:,1));
x = interp1(ro_conv(:,1), ro_conv(:,2), E);
end