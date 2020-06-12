
function In = mom_inertia(type,param)

  In = zeros(3,3);
  switch (type)
    
    % Cylinder: Re, Ri, H
    case 1
      Re = param(1);
      Ri = param(2);
      H  = param(3);
      In(1,1) = ( Re^2 + Ri^2 )/2;
      In(2,2) = ( 3*(Re^2+Ri^2) + H^2 )/12;
      In(3,3) = ( 3*(Re^2+Ri^2) + H^2 )/12;
  
    % Retunr unit matrix
    otherwise
      In = eye(3);
      
  end
  
end

