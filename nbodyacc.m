function [a] = nbodyacc(N, m,r)
%Input arguments
%
%      N:       (real scalar) Number of bodies.
%      m:       (real vector) Length N vector containing masses of bodies.
%      r:       N x 3 array containing the particle positions

%Output:
%      a:        N x 3 array containing the computed particle accelerations

a = zeros(N, 3);

for i = 1:N
    for j = 1:N
        if j ~= i
            res = [0, 0, 0];
            r3ij = ( (r(i,1) - r(j,1))^2 + (r(i,2) - r(j,2))^2 + (r(i,3) - r(j,3))^2)^(3/2);
            res(1) = (m(j) / r3ij) * (r(j,1) - r(i, 1));
            res(2) = (m(j) / r3ij) * (r(j,2) - r(i, 2));
            res(3) = (m(j) / r3ij) * (r(j,3) - r(i, 3));
            a(i,:) = a(i,:) + res;
        end
    end
end
            

end

