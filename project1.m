function project1

delta = calculateDelta(0, 0, 1, 0, 0, 3.2);
fprintf('Feigenbaum constant delta=%d\n', delta);

   function [f, df] = getFunctions(numberOfCycles, r, fInitial, fInitialDerivative)
    % function: getFunctions
    % getFunctions calculates values of derivatives of maps and maps themselves 

    fOld = fInitial;
    fOldDerivative = fInitialDerivative;
    
    iLast = round(2^numberOfCycles);
      
    for i = 1:iLast
      fNew = func(r, fOld);
      fNewDerivative = dFunc(fOld, fOldDerivative);  
      fOld = fNew;
      fOldDerivative = fNewDerivative;
    end 
    
    f = fNew;
    df = fNewDerivative;
    
end;


    function out = calculateRi(r, f0, df0, i)
    % function: calculateRi 
    % calculateRi calculates points of bifurcations
    
    rOld = r;
   
    for n = 1:1000
      [f,df] = getFunctions(i, rOld, f0, df0);
      rNew = rOld - f/df;
      rOld = rNew;
    end 
    
    out = rNew;
    
end;

    function out = calculateDelta(x0, R00, R01, f0, df0, delta)
    % function: calculateDelta 
    % calculateSigma calculates Feigenbaum constant delta for 
    % a discrete map defined by func(r, x)
    % INPUT: 
    % x0 - extremum of a discrete map
    % R00, R01 - Initial approximations for two values of parameter R
    % f0, df0 - Start iterations of function func with f0, and initial 
    %           approximation to the derivative df0 
    % delta - Reasonable initial approximation for the constant
    % OUTPUT:
    % out - value of the constant  
   
    R1 = calculateRi(R00, f0, df0, 0)
    R2 = calculateRi(R01, f0, df0, 1)
    
    for i = 2:10
      R3 = calculateRi(R2+(R2-R1)/delta, f0, df0, i);
      delta = (R2-R1)/(R3-R2);
      R1 = R2;
      R2 = R3;
      fprintf('delta=%d\n', delta);
    end 
    
    out = delta;
    
end;

  
function out = func(r, x)
    % function: func
    % defines a discrete map
    % INPUT: 
    % x -  value of x we want to calculate function at
    % r -  parameter
    % OUTPUT:
    % out - function
         out = r - x^2;
end; % of function func(x)

function out = dFunc(f, df)
    % function: dFunc
    % defines a derivative of a discrete map

    out = 1 - 2*df*f;
end; % of function dFunc(x)


    

end % of project1