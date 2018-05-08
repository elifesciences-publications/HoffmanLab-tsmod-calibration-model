function x = inv_FRETE_dynamic(E, ro)
x = (ro.^6./E - ro.^6).^(1/6);
end
