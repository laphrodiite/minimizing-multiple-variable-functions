function [endVector] = leven_marq(f, fGrad, fHess, xStart, yStart, step, opt)
    k = 1;
    
    grad = matlabFunction(fGrad);
    hess = matlabFunction(fHess);
    numF = matlabFunction(f);
    gradHere = grad(xStart, yStart);
    
    vNew = [xStart yStart];
    I = eye(2);
    
    kMatx = 0;
    xMatx = 0;
    yMatx = 0;
    
    s = linspace(0.01, 2);

    while norm(gradHere) > 0.001
        vOld = vNew;
        x = vOld(1);
        y = vOld(2);
        gradHere = grad(x, y);
        hessianHere = hess(x, y);
        values = eig(hessianHere);
       
        % It is proven that any mk bigger than the bigger eigenvalue of the
        % hessian matrix is acceptable.
        mk = abs(max(values)) + 1;
        dk = linsolve((hessianHere + mk * I), (-1 * gradHere));
        
        if(opt == true)
           x = x + s*dk(1);  y = y + s*dk(2);
           [nim, newStep] = min(numF(x, y));
           step = s(newStep)
        end
        
        vNew(1) = vOld(1) + step * dk(1);
        vNew(2) = vOld(2) + step * dk(2); 
        if numF(vNew(1), vNew(2)) > numF(vOld(1), vOld(2))
            disp("Error: Condition not met.")
        break;
        end

         xMatx(k) = vNew(1);
         yMatx(k) = vNew(2);
         kMatx(k) = k;
        k = k + 1;
    end
    k
    endVector = vNew;

disp(numF(endVector(1), endVector(2)))

    figure(3);
    plot(kMatx, xMatx)
    grid on
    xlabel('k');
    ylabel('xk')
    title('Convergence of x')
    figure(4);
    plot(kMatx, yMatx)
    grid on
    xlabel('k');
    ylabel('yk')
    title('Convergence of y')
end