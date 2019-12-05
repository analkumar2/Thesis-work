function [clr, fitClr] =  getTempColor(tempStr)
if(strcmp(tempStr, '15c'))
    clr = [0 0 1]; fitClr = [1 0 0];
elseif(strcmp(tempStr, '25c'))
    clr = [0 0 0]; fitClr = [1 0 0];
elseif(strcmp(tempStr, '35c'))
    clr = [1 0 0]; fitClr = [0 0 1];
elseif(strcmp(tempStr, 'RT')) % (2013 only)
    clr = [0.1647 0.3843 0.2745]; fitClr = [0 0 1];
else % i.e. 23 nd one 28 (2013 only)
    clr = [0.5843 0.3882 0.3882]; fitClr = [0 0 1];
end