
% Particle Swarm Optimization

function [X,F] = PSO(func,LB,UB,MaxIt,nPop,tol)

  phi1  = 2.05;       % Constriction coefficient 1
  phi2  = 2.05;       % Constriction coefficient 2
  wdamp = 1;          % Inertia Weight Damping Ratio
  
  % Derived parameters
  phi = phi1+phi2;
  chi = 2/(phi-2+sqrt(phi^2-4*phi));
  w = chi;          % Inertia Weight
  c1 = chi*phi1;    % Personal Learning Coefficient
  c2 = chi*phi2;    % Global Learning Coefficient

  % Limits/Integer variables
  nVar = length(LB);
  pLB = repmat(LB,nPop,1);
  pUB = repmat(UB,nPop,1);
  VelMax = 0.1*(pUB-pLB);
  VelMin = -VelMax;

  % Initialize
  GlobalBestCost = inf;
  particlePosition = unifrnd(pLB,pUB);
  particleVelocity = zeros(nPop,nVar);
  particleCost = zeros(nPop,1);
  for i = 1:nPop  
    particleCost(i) = func( ( particlePosition(i,:) )' );
  end
  particleBestPosition = particlePosition;
  particleBestCost = particleCost;

  % Update Global Best
  for i = 1:nPop
    if ( particleBestCost(i) <= GlobalBestCost )
      GlobalBestCost = particleBestCost(i);
      GlobalBestPosition = particleBestPosition(i,:);
    end    
  end

  % Main Loop
  Ftmp = zeros(MaxIt,1);
  for it = 1:MaxIt

    % Update Velocity
    particleVelocity = w*particleVelocity ...
          +c1*rand(nPop,nVar).*(particleBestPosition-particlePosition) ...
          +c2*rand(nPop,nVar).*(repmat(GlobalBestPosition,nPop,1)-particlePosition);
        
    % Apply Velocity Limits
    particleVelocity = max(particleVelocity,VelMin);
    particleVelocity = min(particleVelocity,VelMax);
        
    % Update Position
    particlePosition = particlePosition + particleVelocity;
    
    % Velocity Mirror Effect
    IsOutside = (particlePosition<pLB | particlePosition>pUB);
    particleVelocity(IsOutside) = -particleVelocity(IsOutside);
        
    % Apply Position Limits
    particlePosition = max(particlePosition,pLB);
    particlePosition = min(particlePosition,pUB);
        
    % Evaluation
    for i = 1:nPop  
      particleCost(i) = func( ( particlePosition(i,:) )' );
    end
        
    % Update Personal Best
    flag = (particleCost < particleBestCost);       
    particleBestCost(flag) = particleCost(flag);
    particleBestPosition(flag,:) = particlePosition(flag,:);

    % Update Global Best
    for i = 1:nPop  
      if ( particleBestCost(i) <= GlobalBestCost )
        GlobalBestCost = particleBestCost(i);
        GlobalBestPosition = particleBestPosition(i,:);
      end  
    end
    
    w = w*wdamp;
    
    % Save intermediate cost function
    Ftmp(it,1) = GlobalBestCost;
    
    if ( GlobalBestCost < tol )
      break;
    end
    
  end

  X = GlobalBestPosition;
  F = Ftmp(1:it,1);

end     % End of function

