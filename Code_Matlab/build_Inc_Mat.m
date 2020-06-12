
% Function to build the incidence matrices SS, S0, SE
% See manual page 9
% by GabGil 23 July 2004

function [S0,SS,SE] = build_Inc_Mat(BB)

% Check vector BB and its elements
n = length(BB);
if (n == 0)
  fprintf('\n No elements in matrix BB\n');
  return
end
for j = 1:n
  if ((BB(j) < 0) || (BB(j) > n)) 
    fprintf('\n Element BB(%d) out of bounds\n',j);
    return
  end
end

% Build vector S0
S0 = zeros(1,n);
for j = 1:n
  idx = BB(j);
  if (idx == 0)
    S0(1,j) = 1;
  else
    S0(1,j) = 0;
  end
end

% Build matrix SS
SS = zeros(n,n);
for j = 1:n
  idx = BB(j);
  for i = 1:n
    if (idx == i) 
      SS(i,j) = 1;
    elseif (i == j) 
      SS(i,j) = -1;
    else
      SS(i,j) = 0;
    end
  end
end

% Build vector SE
SE = zeros(1,n);
for i = 1:n
  flag = 0;
  for j = 1:n
    if (SS(i,j) == 1)
      flag = 1;
    end
  end
  if (flag == 0)
    SE(1,i) = 1;
  else
    SE(1,i) = 0;
  end
end

% EOF
