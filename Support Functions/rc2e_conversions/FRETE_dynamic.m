%% Function to Convert End to End Distance to FRET Efficiency
function E = FRETE_dynamic(x, ro)
E = ro.^6./(ro.^6+x.^6);
end
