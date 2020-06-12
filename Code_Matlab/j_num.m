%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	May 3, 1998, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   J_NUM	Find joint connection from a given End-link to the 0-th link
%
%		1998.1.9  A.Kurosu
%

function joint = j_num( num_e )

global BB SE

n = length(SE);

j = 0;

for i = 1 : n
   
   if SE(i) == 1
      j = j + 1;
      ie(j) = i;
      
   end
   
end

j_number = BB(ie(num_e));
connection = [ie(num_e)];

while (j_number ~= 0)
   
   connection = [j_number connection];
   j_number = BB(j_number);
   
end

joint = connection;


%%%EOF
