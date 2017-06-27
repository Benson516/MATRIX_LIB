
clc

%%
% Generate the M.txt file
 
% M = [1 2 0;
%      4 5 6;
%      0 8 9]
 
%  M = [1 2 3 4 5;
%      6 7 8 9 0;
%      5 4 3 2 1;
%      0 9 8 7 6;
%      -1 1 -1 1 -1]

% M = [1 0 2 0;
%      0 1 0 0;
%      0 0 0 1;
%      0 0 1 0]

M = rand(100,100);

rank_M = rank(M)

is_invertable = rank_M == size(M,1)

M_inv = M^-1;

if size(M,1) <= 10
    M_inv
end

writeMatrixToFile('./M.txt', M);
