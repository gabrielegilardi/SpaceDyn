
function [ PF0, PFe, Ptau ] = calc_Work( v0, w0, vve, ww, qd, F0, T0, Fe, ...
                                    Te, tau )

  num_q = length(qd);
  
  %Work forces/moments on the base
  PF0 = F0'*v0 + T0'*w0;
  
  %Work forces/moments on the end effector
  PFe = zeros(1,num_q);
  for i = 1:num_q
    PFe(1,i) = Fe(:,i)'*vve(:,i) + Te(:,i)'*ww(:,i);
  end
  
  %Work forces/moments on the joints
  Ptau = zeros(1,num_q);
  for i = 1:num_q
    Ptau(1,i) = tau(i)*qd(i);
  end

end     % End of function
