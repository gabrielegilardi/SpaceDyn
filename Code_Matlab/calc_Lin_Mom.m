
function [ LM0, LMj ] = calc_Lin_Mom( v0, vv )

  global m0 m

  num_q = size(vv,2);
  LM0  = m0*v0;
  LMj  = zeros(3,num_q);
  for i = 1:num_q
    LMj(:,i)  = m(i)*vv(:,i);
  end

end     % End of function
