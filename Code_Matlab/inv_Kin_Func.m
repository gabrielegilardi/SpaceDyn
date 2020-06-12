
function  f = inv_Kin_Func(x)

  global paramPSO

  POS = paramPSO.POS;
  type = paramPSO.type;
  
  % Full variables
  R0 = x(1:3);
  Q0 = x(4:6);
  q  = x(7:end);

  % Common quantities
  A0 = rpy2dc(Q0)';
  AA = calc_aa( A0, q );
  RR = calc_pos( R0, A0, AA, q );

  % IK on the center of mass
  if ( type == 0 )
    now_p = calc_pos_CoM( R0, RR);
  
  % IK on an end-effector
  else
    joints = j_num( type );
    [now_p, ~] = f_kin_e( RR, AA, joints );
  end

  % Error function
  err = POS - now_p;
  f = sum( err.^2 )/2;

end
