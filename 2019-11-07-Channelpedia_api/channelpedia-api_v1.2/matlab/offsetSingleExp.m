function [xIn, yOut, strEq] = offsetSingleExp(xIn, C, A, Tau, expType)
% singleExp returns single exponential value for a given A, Tau and expType
% type expType
% Usage   : [xIn, yOut, strEq] = singleExp(xIn, yIn, Xout, 1)  
% expType : 1 for C+( A * exp(-t/Tau))   ,  Left  Down/concave
%         : 2 for C+( A * exp( t/Tau))   ,  Right Down/concave
%         : 3 for C+(A - exp( t/Tau))   ,  Left  Up  /Convex
%         : 4 for C+(A *(1.0-exp(-t/Tau))),Right Up  /Convex

if(Tau == 0)
    Tau = 1e-99;
end
switch expType
    case 1
        yOut = C + (A .* exp(-xIn ./ Tau));
        strEq = sprintf('%.3f + (%.3f * exp(-t/%.3f))', C, A, Tau);
    case 2
        yOut = C + (A .* exp(xIn ./ Tau));
        strEq = sprintf('%.3f + (%.3f * exp(t/%.3f))', C, A, Tau);
    case 3
        yOut = C + (A  - exp(xIn ./ Tau));
        strEq = sprintf('%.3f + (%.3f - exp(t/%.3f))', C, A, Tau);
    case 4
        yOut = C + (A .*(1.0 - exp(-xIn ./ Tau)));
        strEq = sprintf('%.3f + (%.3f *(1.0 - exp(-t/%.3f)))', C, A, Tau);
    otherwise
        yOut = [];
        xIn  = [];
        strEq = 'Invalide expType';
end