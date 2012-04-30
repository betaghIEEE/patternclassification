% Test the forward dynamic programming algorithm
% Copyright 1999 by Todd K. Moon

G{1} = [2,3,4,5];  G{2} = [6,7];  G{3} = [7,8];  G{4} = [6,8];
G{5} = 8;   G{6} = [9,10,11];  G{7} = [10,11,12];  G{8} = 11;
G{9} = 13;  G{10} = 13;  G{11} = 13;  G{12} = 13;
W{1} = [6,2,3,4];  W{2} = [2,3];  W{3} = [9,4];  W{4} = [3,5];
W{5} = 2;   W{6} = [2,3,7];  W{7} = [5,2,6];  W{8} = 3;
W{9} = 4;   W{10} = 2;   W{11} = 3;   W{12} = 5;

[plist,cost] = fordyn(G,W)