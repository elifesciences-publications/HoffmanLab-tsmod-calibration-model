%% function to convert P(re) to P(rc)
function rc = re_2_rc(x, conv)
if max(x) > max(conv(:,2)) 
    disp('Need longer vector in conv')
end
rc = interp1(conv(:,1), conv(:,2), x);
end