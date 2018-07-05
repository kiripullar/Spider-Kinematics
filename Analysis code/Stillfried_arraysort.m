function L = Stillfried_arraysort(M,varargin)
% INDEX creates a list with indices and values of an array.
%
% SYNTAX:
% L = index(M);
% L = index(M,'option');
% L = index(M,'option1','option2',...);
% 
% DESCRIPTION:
% The function index.m takes a matrix or an array and creates a list that 
% contains the indices and the values of the matrix. 
% In the resulting matrix L, there is one row for each cell of M.
% The first elements of a row are the index, the last one is the value of the cell.
% 
% OPTIONS:
% There are three possible options: 'sorted', 'nonans' and 'nozeros'. 
% The option 'sorted' returns a list that is sorted according to the
% values.
% The option 'nonans' returns a list with NaNs ignored.
% With the option 'nozeros' zeros are ignored.
%
% The options can be entered in any order.
%
% Example:
%
% >> a = [2,3,NaN;1,2,4;7,0,9]
%
%      a =
%      2 3 NaN
%      1 2 4
%      7 0 9
%
% >> i = index (a,'sorted','nozeros')
%
%     i =
%     2 1 1
%     1 1 2
%     2 2 2
%     1 2 3
%     2 3 4
%     3 1 7
%     3 3 9
%     1 3 NaN
% 
% 
% REMARK:
% When non-integer numbers are used, the indices are displayed as floats, but stored as integers:
% 
% >> a=[2.2,3.1,NaN;1.9,2.3,4;7,0,9]
% 
% a =
% 
%     2.2000    3.1000       NaN
%     1.9000    2.3000    4.0000
%     7.0000         0    9.0000
% 
% >> i = index(a,'sorted','nonans','nozeros')
% 
% i =
% 
%     2.0000    1.0000    1.9000
%     1.0000    1.0000    2.2000
%     2.0000    2.0000    2.3000
%     1.0000    2.0000    3.1000
%     2.0000    3.0000    4.0000
%     3.0000    1.0000    7.0000
%     3.0000    3.0000    9.0000
% 
% >> i(:,1:2)
% 
% ans =
% 
%      2     1
%      1     1
%      2     2
%      1     2
%      2     3
%      3     1
%      3     3
%
% version 1.3: speed improvement by using ind2sub
% history:
%  version 1.2: speed improvement by using repmat
%  version 1.1: some code simplification, some change in the help text
% 
% (C) 2008-2009 Georg Stillfried (gsti@gmx.net)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S=size(M);
dims=ndims(M);

% make list
[ind{1:dims}] = ind2sub(S,1:numel(M));  
L = [cat(1,ind{:}).' M(:)];  

% options
for i = 1:nargin-1
    var = varargin{i};
    if ischar(var)
        switch var
            case 'nonans'
                %remove lines where value == NaN
                L = L(~isnan(L(:,dims+1)),:);
            case 'nozeros'
                %remove lines where value == 0
                L = L(L(:,dims+1)~=0,:);
            case 'sorted'
                %sort
                [b,ix] = sort(L(:,dims+1));
                L = L(ix,:);
        end
    end
end