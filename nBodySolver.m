function [t, pos] = nBodySolver(tmax, level, N, masses, pos0, v0, tracefreq)

%  Solves and plots nbody problem using FDA 
%  
%  Input arguments
%
%      tmax:       (real scalar) Final solution time.
%      level:      (integer scalar) Discretization level.
%      N:          (real scalar) Number of bodies.
%      masses:     (real vector) Length N vector containing masses of bodies.
%      pos0:       (real matrix) Nx3 matrix with initial position of each
%                                body
%      v0:         (real matrix) Nx3 matrix with initial velocity of each
%                                body
%      tracefreq:  (optional integer scalar) Frequency of tracing output, 
%                  0 disables tracing.
%
%  Output arguments
%
%      t:      (real vector) Vector of length nt = 2^level + 1 containing
%              discrete times (time mesh).
%      pos:    N x 3 x nt array giving position of every body at each
%              discrete time
%     

% Tracing control: if 5th arg is supplied base tracing on that input,
% otherwise use local defaults.
   if nargin > 6
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
   
   % Define number of time steps and create t and pos arrays
   nt = 2^level + 1;
   t = linspace(0.0, tmax, nt);
   pos = zeros(N, 3, nt);
   
    % Determine discrete time step from t array.
   deltat = t(2) - t(1);
   
   % Initialize first two values of the each body's position
   pos(:,:,1) = pos0;
   pos(:,:,2) = pos0 + deltat*v0 + 0.5*deltat^2 * nbodyacc(N, masses, pos0);
   
   % Evolve the body position through every time step from n=3 to n=nt
   for n = 2 : nt-1
       if rem(n, tracefreq) == 0
         fprintf('nbody: Step %g of %g\n', n, nt);
       end
       
       pos(:,:,n+1) = 2*pos(:,:,n) - pos(:,:,n-1) + (deltat^2) *  nbodyacc(N, masses, pos(:,:,n));

   end

   
% Create an animated plot of the xy coordinates of both masses over time
plotenable = 1;
pausesecs = 0.0;
n = 1;
dlim = 10;
% for ti = 0: deltat : tmax
%     
%     if plotenable
%         clf;
%         hold on;
%         axis square;
%         box on;
%         
%         xlim([-dlim, 1 + dlim]);
%         ylim([-dlim, 1 + dlim]);
%         plot(pos(1,1,n), pos(1,2,n), 'Marker', 'o', 'MarkerSize', 15, ...
%          'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
%         plot(pos(2,1,n), pos(2,2,n), 'Marker', 'o', 'MarkerSize', 15, ...
%          'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
% 
% 
%         
%         drawnow;
%         pause(pausesecs);
%     end
%     n = n+1;
% end

end
   
   

   
   
   
   





