
function [ HM0d, HMjd ] = calc_der_Ang_Mom( R0, A0, RR, AA, v0, w0, ...
                                            vv, ww, vd0, wd0, vd, wd, ...
                                            Pref, Vref)

  global m0 m inertia0 inertia

  num_q = size(vv,2);
  In = A0*inertia0*A0';
  HM0d = inertia0*wd0 + cross( w0, In*w0 ) + cross( v0-Vref, m0*v0 ) + ...
         cross( R0-Pref, m0*vd0 );
  HMjd = zeros(3,num_q);
  for i = 1:num_q
    In = AA(:,i*3-2:i*3)*inertia(:,i*3-2:i*3)*AA(:,i*3-2:i*3)';
    HMjd(:,i) = In*wd(:,i) + cross( ww(:,i), In*ww(:,i) ) ...
                + cross( vv(:,i)-Vref, m(i)*vv(:,i) ) ...
                + cross( RR(:,i)-Pref, m(i)*vd(:,i) ); 
  end
  
end   % End of function
