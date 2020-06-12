
function [ HM0, HMj ] = calc_Ang_Mom( R0, A0, RR, AA, v0, w0, vv , ww, Pref )

  global m0 m inertia0 inertia

  num_q = size(vv,2);
  In = A0*inertia0*A0';
  HM0  = cross( R0-Pref , m0*v0 ) + In*w0;
  HMj  = zeros(3,num_q);
  for i = 1:num_q
    In = AA(:,i*3-2:i*3)*inertia(:,i*3-2:i*3)*AA(:,i*3-2:i*3)';
    HMj(:,i)  = cross( RR(:,i)-Pref , m(i)*vv(:,i)) + In*ww(:,i);
  end
  
end   % End of function
