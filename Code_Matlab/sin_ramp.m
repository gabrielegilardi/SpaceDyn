
function  f = sin_ramp(t,dt)

 y = t/dt - sin(2*pi*t/dt)/(2*pi); 
 
 f = min(y,1);

end

