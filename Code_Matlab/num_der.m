
function der = num_der(f,t)

  n = length(f);
  der = zeros(n,1);
  if ( length(t) == 1 )
      dt = t;
      t = 0:dt:(n-1)*dt;
  end
  % 1st element always zero
  der(1) = 0;
  if ( n == 1 )
    return
  end
  % 2nd element using 1st order left-sided equation
  der(2) = ( f(2) - f(1) )/(t(2)-t(1));
  if ( n == 2 )
    return
  end
  % All other elements using a 2nd order left-sided equation
  for i = 3:n
    der(i) = ( 3*f(i) - 4*f(i-1) + f(i-2) ) / (t(i) - t(i-2));
  end

end
