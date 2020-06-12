
function CoMv = calc_vel_CoM( v0, vv)

global m0 m

  CoMv = (m0*v0 + sum(vv*diag(m),2))/(m0 + sum(m));

end

