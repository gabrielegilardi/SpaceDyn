%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1.2 // Sep. 11, 2002, Last modification by H.Hamano
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% CALC_JE	Calculation of the Jacobian Matrix (6 by n)
%           for the endpoint given by connection vector 'joints.'
%
%     1997.1.14 Space Robotics Lab (C)
%               A.Kurosu
%     1998.2.4  K.Yoshida
%     2002.9.11 H.Hamano
%

function Jacobian = calc_je( RR, AA, q, joints )

num_q = length(q);

% number of links from base to endpoint.
n = length(joints);

% Calculation of Jacobian
% JJ_te = zeros(3,num_q);
% JJ_re = zeros(3,num_q);

JJ_te = calc_jte( RR, AA, q, joints );
JJ_re = calc_jre( AA, joints );
JJ = [ JJ_te; JJ_re ];

% Compose the Jacobian using the corresponding joints.
Jacobian = zeros(6,num_q);

for i = 1 : 1 : n
   Jacobian(1:6 , joints(i)) = JJ(1:6 , i);
end


%%%EOF
