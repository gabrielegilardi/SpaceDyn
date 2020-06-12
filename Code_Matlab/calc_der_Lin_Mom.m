
function [ LM0d, LMjd ] = calc_der_Lin_Mom( vd0, vd )

  global m0 m

  num_q = size(vd,2);
  LM0d = m0*vd0;
  LMjd = zeros(3,num_q);
  for i = 1:num_q
    LMjd(:,i) = m(i)*vd(:,i);
  end

end     % End of function
