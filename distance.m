function dist = distance(parameters, X, Y)
    size_X = size(X);
    size_X = size_X(1);
    dist = NaN(size_X,1);
    for i = 1:size_X
        dists = zeros(4,1);
        dists(1) = abs(parameters(1) - X(i));
        dists(2) = abs(parameters(2) - Y(i));
        dists(3) = abs(parameters(3) - X(i));
        dists(4) = abs(parameters(4) - Y(i));
        dist(i) = min(dists);
    end
return
