%% function to convert P(re) to P(rc)
function re = rc_2_re(x, conv)
if max(x) > max(conv(:,2)) 
    disp('Need longer vector in conv')
end
re = interp1(conv(:,2), conv(:,1), x);
end