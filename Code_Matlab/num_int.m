
function int = num_int(f,dt)

  n = length(f);
  int = zeros(n,1);
  for i = 2:n
    int(i) = int(i-1) + dt*( f(i-1) + f(i) )/2;
  end

end

