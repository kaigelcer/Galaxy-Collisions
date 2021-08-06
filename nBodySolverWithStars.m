function [outputArg1,outputArg2] = nBodySolverWithStars(tmax, level, masses, pos0, v0, nStars, rotationDirection, avienable, tracefreq)

% "Constants"
X = 1;
Y = 2;
Z = 3;
N = 2; % Number of cores

%  Simulate Galaxy Collisions by building on the nbody solver
%  Input arguments
%
%      tmax:       (real scalar) Final solution time.
%      level:      (integer scalar) Discretization level.
%      masses:     (real vector) Length 2 column vector containing masses of core.
%      pos0:       (real matrix) 2x3 matrix with initial position of each
%                                core
%      v0:         (real matrix) 2x3 matrix with initial velocity of each
%                                core
%      nStars:     number of total massless stars, to be split evenly
%                  amongst cores
%      rotationDirection: enter 1 for clockwise rotation of stars, 2 for
%      counterclockwise
%      avienable:  1 to enable avi creation, 0 to disable
%      tracefreq:  (optional integer scalar) Frequency of tracing output, 
%                  0 disables tracing.
%
%  Output arguments
%
%      t:      (real vector) Vector of length nt = 2^level + 1 containing
%              discrete times (time mesh).
%      pos:    2 x 3 x nt array giving position of each core at each
%              discrete time
%     

% Tracing control: if 5th arg is supplied base tracing on that input,
% otherwise use local defaults.
   if nargin > 8
      if tracefreq == 0
         trace = 0;
      else
         trace = 1;
      end
   else
      trace = 1;
      tracefreq = 100;
   end

   if trace
      fprintf('In nbody: Argument dump follows\n');
      tmax, level, masses, pos0, v0
   end
   
   
   % Create matrix of star initial positions, randomly distributing them
   % about each core
   starsPos0 = zeros(nStars, 3);
   maxRadius = 0.5;
   for star = 1 : 1 : nStars/2
       starsPos0(star, X) = 2*maxRadius*rand - maxRadius + pos0(1,X);
       starsPos0(star, Y) = rand*2*sqrt(maxRadius^2 - (starsPos0(star,X) - pos0(1,X))^2)-sqrt(maxRadius^2 - (starsPos0(star,X) - pos0(1,X))^2) + pos0(1,Y);
   end
   for star = nStars/2 + 1 : 1 : nStars
       starsPos0(star, X) = 2*maxRadius*rand - maxRadius + pos0(2,X);
       starsPos0(star, Y) = rand*2*sqrt(maxRadius^2 - (starsPos0(star,X) - pos0(2,X))^2)-sqrt(maxRadius^2 - (starsPos0(star,X) - pos0(2,X))^2) + pos0(2,Y);
   end
   
   %Set star initial velocities to that of the core that they're orbiting
   %Also create a radii vector of the initial distance of each star from
   %the core
   starsv0 = zeros(nStars, 3);
   starsRadii = zeros(nStars, 3);
   for star = 1 : 1 : nStars/2
       starsv0(star, X) = v0(1,X);
       starsv0(star, Y) = v0(1,Y);
       starsRadii(star, X) = starsPos0(star, X) - pos0(1,X);
       starsRadii(star, Y) = starsPos0(star, Y) - pos0(1,Y);
   end
   
   for star = nStars/2 + 1 : 1 : nStars
       starsv0(star, X) = v0(2,X);
       starsv0(star, Y) = v0(2,Y);
       starsRadii(star, X) = starsPos0(star, X) - pos0(2,X);
       starsRadii(star, Y) = starsPos0(star, Y) - pos0(2,Y);
   end
      
   %Setting initial velocity such that stars orbit
   for star = 1:1:nStars
       pos0Temp = [pos0; starsPos0(star,:)]; %Creates an array with both core positions and 1 star's position
       massesTemp = [masses;0];
       acceleration = nbodyacc(N+1,massesTemp,pos0Temp); % Calculate acceleration
       starAcceleration = acceleration(3,:);
       
       %Get velocity direction and unit vector for orbit
       if rotationDirection == 1
           velocityDirection = [starsRadii(star, Y), - starsRadii(star, X), 0];
       else
           velocityDirection = [- starsRadii(star, Y), starsRadii(star, X), 0];
       end
       velocityUnitVector = velocityDirection / norm(velocityDirection);
       
       %Calculates the velocity required for orbit
       starsv0(star, :) = starsv0(star,:) + (velocityUnitVector) * sqrt(norm(starAcceleration)*norm(starsRadii(star,:)));
       
   end
   
   
   % Define number of time steps and create t and pos arrays
   nt = 2^level + 1;
   t = linspace(0.0, tmax, nt);
   pos = zeros(N, 3, nt);
   starsPos = zeros(nStars, 3, nt);
   
    % Determine discrete time step from t array.
   deltat = t(2) - t(1);
   
   % Initialize first two values of the each core's position
   pos(:,:,1) = pos0;
   
   pos(:,:,2) = pos0 + deltat*v0 + 0.5*deltat^2 * nbodyacc(N, masses, pos0);
   
   % Initialize first two values of the each star's position
   starsPos(:,:,1) = starsPos0;
   for star = 1:1:nStars
       pos0Temp = [pos0; starsPos0(star,:)]; %Creates an array with both core positions and 1 star's position
       v0Temp = [v0; starsv0(star,:)]; %Creates an array with both core velocities and 1 star's velocity
       massesTemp = [masses; 0]; %Adds an additional 0 mass to masses vector
       posTemp = pos0Temp + deltat*v0Temp + 0.5*deltat^2 * nbodyacc(N+1, massesTemp, pos0Temp); 
       starsPos(star,:,2) = posTemp(3,:);
   end
      
   for n = 2 : nt-1
       if rem(n, tracefreq) == 0
         fprintf('nbody: Step %g of %g\n', n, nt);
       end
       
       %Calculate next position of every star
       for star = 1:1:nStars
           posnTemp = [pos(:,:,n);starsPos(star,:,n)]; %Creates an array with both core positions and 1 star's position
           posnMinus1Temp = [pos(:,:,n-1);starsPos(star,:,n-1)]; %Same as above, but for step n-1
           massesTemp = [masses; 0];
           %Calculate star's next position
           posnPlus1Temp = 2*posnTemp - posnMinus1Temp + (deltat^2) * nbodyacc(N+1, massesTemp, posnTemp);
           starsPos(star,:,n+1) = posnPlus1Temp(3,:);
       
       %Calculate next position of each core
       pos(:,:,n+1) = 2*pos(:,:,n) - pos(:,:,n-1) + (deltat^2) *  nbodyacc(N, masses, pos(:,:,n));
       end
       
   end


% ANIMATE AND/OR PLOT

%For movie creation
avifilename='galaxy.avi';
if avienable
   aviobj = VideoWriter(avifilename);
   open(aviobj);
end
aviframerate = 25;


plotenable = 1;
pausesecs = 0;
n = 1;
dlim = 10;
for t = 0: deltat : tmax
    
    if plotenable
        clf;
        hold on;
        axis square;
        box on;
        
        xlim([-dlim, 1 + dlim]);
        ylim([-dlim, 1 + dlim]);
        
        plot(pos(1,X,n), pos(1,Y,n), 'Marker', 'o', 'MarkerSize', 3, ...
        'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
        
        plot(pos(2,X,n), pos(2,Y,n), 'Marker', 'o', 'MarkerSize', 3, ...
        'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
        
    
        scatter(starsPos(1:nStars/2,X,n), starsPos(1:nStars/2,Y,n), 1, 'b', 'o', 'Filled');
    
        scatter(starsPos(nStars/2:nStars,X,n), starsPos(nStars/2:nStars,Y,n), 1, 'g', 'o', 'Filled');
    
        
        drawnow;
        pause(pausesecs);
    end
    
    if t == 0
        framecount = 5 * aviframerate ;
    else
        framecount = 1;
    end
    for iframe = 1 : framecount
        writeVideo(aviobj, getframe(gcf));
    end
         
    
    n = n+1;
end

end
