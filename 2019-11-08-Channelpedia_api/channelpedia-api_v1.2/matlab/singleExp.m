function [xIn, yOut, strEq] = singleExp(xIn, A, Tau, expType)
% singleExp returns single exponential value for a given A, Tau and expType
% type expType
% Usage   : [xIn, yOut, strEq] = singleExp(xIn, yIn, Xout, 1)  
% expType : 1 for A * exp(-t/Tau)   ,  Left  Down/concave
%         : 2 for A * exp( t/Tau)   ,  Right Down/concave
%         : 3 for A - exp( t/Tau)   ,  Left  Up  /Convex
%         : 4 for A *(1.0-exp(-t/Tau)),Right Up  /Convex

if(Tau == 0)
    Tau = 1e-99;
end
switch expType
    case 1
        yOut = A .* exp(-xIn ./ Tau);
        strEq = sprintf('%.3f * exp(-t/%.3f)', A, Tau);
    case 2
        yOut = A .* exp(xIn ./ Tau);
        strEq = sprintf('%.3f * exp(t/%.3f)', A, Tau);
    case 3
        yOut = A  - exp(xIn ./ Tau);
        strEq = sprintf('%.3f - exp(t/%.3f)', A, Tau);
    case 4
        yOut = A .*(1.0 - exp(-xIn ./ Tau));
        strEq = sprintf('%.3f *(1.0 - exp(-t/%.3f))', A, Tau);
    otherwise
        yOut = [];
        xIn  = [];
        strEq = 'Invalide expType';
end