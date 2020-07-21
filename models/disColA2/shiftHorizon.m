function u0 = shiftHorizon(u)
%SHIFTHORIZON Summary of this function goes here
% 
% [OUTPUTARGS] = SHIFTHORIZON(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2019/01/15 15:40:58 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2019

u0 = [u(:,2:size(u,2)) u(:,size(u,2))];

end
