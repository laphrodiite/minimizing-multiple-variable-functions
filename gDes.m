function [endVector] = gDes(f, fGrad, xStart, yStart, step, opt)

syms x y
k = 1;

numF = matlabFunction(f);
grad = matlabFunction(fGrad);

vNew(1) = xStart
vNew(2) = yStart
x = vNew(1);
y = vNew(2);
gradHere = grad(x, y);

kMatx = 0;
xMatx = 0;
yMatx = 0;

s = linspace(0.01, 2);


while norm(gradHere)>0.001
    vOld = vNew;
    x = vOld(1);
    y = vOld(2);
    gradHere = grad(x, y);
    dk = (-1) * gradHere;
        
    % For the second part
    if(opt == true)
       x = x + s*dk(1);  y = y + s*dk(2);
       [nim, newStep] = min(numF(x, y));
       step = s(newStep)
    end
       
    vNew(1) = vOld(1) + step*dk(1);
    vNew(2) = vOld(2) + step*dk(2);
    
    if numF(vNew(1), vNew(2)) > numF(vOld(1), vOld(2))
        disp("Error: Condition not met.")
        break;
    end

    xMatx(k) = vNew(1);
    yMatx(k) = vNew(2);
    kMatx(k) = k;
    k = k + 1;
end

endVector = vNew
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
