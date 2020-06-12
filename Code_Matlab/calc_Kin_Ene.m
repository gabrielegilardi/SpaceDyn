
function [ TK0, TKj ] = calc_Kin_Ene( A0, AA, v0, w0, vv, ww )

  global m0 m inertia inertia0

  num_q = size(vv,2);
  TK0 = m0*( v0'*v0 )/2 + w0'*A0*inertia0*A0'*w0/2;
  TKj = zeros(1,num_q);
  for i = 1:num_q
    An = AA(:,i*3-2:i*3);
    In = inertia(:,i*3-2:i*3);
    TKj(i) = m(i)*( vv(:,i)'*vv(:,i) )/2 + ww(:,i)'*An*In*An'*ww(:,i)/2;
  end

end     % End of function

