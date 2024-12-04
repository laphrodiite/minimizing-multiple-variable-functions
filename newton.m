function [endVector] = newton(f, fGrad, fHess, xStart, yStart, step, minf, armijo)
    
    syms x y;
    k = 1;
    
    grad = matlabFunction(fGrad);
    numF = matlabFunction(f);
    hess = matlabFunction(fHess);
    
    vNew(1) = xStart;
    vNew(2) = yStart;
    x = vNew(1);
    y = vNew(2); 
    gradHere = grad(x, y);
    
    kMatx = 0;
    xMatx = 0;
    yMatx = 0;
    
    s = linspace(0.01, 2);
    
    % Armijo parameters
    alpha = 0.5; % Reduction factor
    beta = 0.1;  % Armijo threshold
    
    while norm(gradHere)>0.001 & k < 1000
        vOld = vNew;
        x = vOld(1);
        y = vOld(2);
        gradHere = grad(x, y);
        hessianHere = hess(x, y);
        dk = hessianHere\gradHere; %inv(hessianHere) * gradHere
       
        % For minimizing f
        if(minf == true)
               x = x + s*dk(1);  y = y + s*dk(2);
               [nim, newStep] = min(numF(x, y));
               step = s(newStep);
        end
    
        % For Armijo rule
        if (armijo == true)
            step = 1; % Start with an initial step size
            iter = 0;
            while (numF(x + step * dk(1), y + step * dk(2)) > ...
                  numF(x, y) + beta * step * (gradHere * dk') & iter < 500)
                step = alpha * step; % Reduce step size
                iter = iter + 1;
            end
        end
       
        vNew(1) = vOld(1) - step * dk(1);
        vNew(2) = vOld(2) - step * dk(2);
        
        if numF(vNew(1), vNew(2)) > numF(vOld(1), vOld(2))
            disp("Error: Condition not met.")
        end
        xMatx(k) = vNew(1);
        yMatx(k) = vNew(2);
        kMatx(k) = k;
    
        k = k +1;
    end
    
    disp("Newton:")
    endVector = vNew
    disp(numF(vNew(1), vNew(2)))
    
    figure(5); % Single figure for both plots
    
    % Plot 1: Convergence of x
    subplot(2, 1, 1); % First subplot (2 rows, 1 column, 1st position)
    plot(kMatx, xMatx);
    grid on;
    xlabel('k');
    ylabel('x_k');
    title('Convergence of x in [1,-1]');
    
    % Plot 2: Convergence of y
    subplot(2, 1, 2); % Second subplot (2 rows, 1 column, 2nd position)
    plot(kMatx, yMatx);
    grid on;
    xlabel('k');
    ylabel('y_k');
    title('Convergence of y in [1,-1]');
    sgtitle('Newton');

end