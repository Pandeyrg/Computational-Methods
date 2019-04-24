function tc = orthoProjectionOnCurve(x0, y0, X, Y, dXdt, dYdt, eps)
    fun =@(t) dXdt(t)*( X(t)-x0 ) + dYdt(t)*(Y(t) - y0) ; %first derivative of the distance function .
    dist = @(t) (X(t)- x0)*(X(t) -x0) + ( Y(t) - y0)*(Y(t) - y0); % distance function
    % we call the point (x0, y0 ) as the reference point.
    tc =0;
    min_t = 0; % to store the root of f'(t) corresponding to minimum distance so far.
    min_dis =inf; % to store the minimum distance of the reference point from all the roots of f'(t) obtained so far.
    max_trials = 500; % maximum no of trials for Secant Method.
    part = 1024; % gives the number of parts in which we segment the curve.
    
    for k = 1: part   %iterating over all the parts of the curve.
        t1 = k/part;
        t2 = (k+1)/part;
        if t2 > 1 
            t2 = t2 - 1;
        end
        f_1 =  fun(t1);
        trials = 0;
        d = inf;
        err = inf;
        f_1 = fun(t1);
        while ( err > eps && trials < max_trials ) % Secant Algorithm
            f_2 = fun(t2);
            d = (f_2*(t2 - t1))/(f_2 - f_1);
            t1 =  t2;
            t2 =  t2 - d;
            f_1  = f_2;
            err =  dXdt(t2)*(X(t2) - x0) + dYdt(t2)*(Y(t2) -y0); 
            err =  err / sqrt( (dXdt(t2)*dXdt(t2)+dYdt(t2)*dYdt(t2))*((X(t2) - x0)*(X(t2) - x0)+(Y(t2) -y0)*(Y(t2) -y0)));
            % err is the cosine of the angle betweent the tangent vector at
            % t2 on the curve and the projection vector
        end
        dis = dist(t2);
         if dis < min_dis % checking if we find a root with a lesser distance than computed earlier.
             min_t = t2;
             min_dis = dis;
         end
    end
    tc = min_t-floor(min_t); % returns tc with value in between 0 and 1
end
    

   