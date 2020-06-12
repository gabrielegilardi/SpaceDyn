
function CoM = calc_pos_CoM( R0, RR)

global m0 m

  CoM  = (m0*R0 + sum(RR*diag(m),2))/(m0 + sum(m));

end

