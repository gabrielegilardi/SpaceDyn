
function CoMa = calc_acc_CoM( vd0, vd)

global m0 m

  CoMa = (m0*vd0 + sum(vd*diag(m),2))/(m0 + sum(m));

end

