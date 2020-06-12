
function [ VG0, VGj ] = calc_Pot_Ene( R0, RR )

  global m0 m Gravity

  num_q = size(RR,2);
  VG0 = -m0*( Gravity'*R0 );
  VGj = zeros(1,num_q);
  for i = 1:num_q
    VGj(i) = -m(i)*( Gravity'*RR(:,i) );
  end

end     % End of function
