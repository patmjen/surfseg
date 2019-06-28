function S = RobustSign(X,tol)
    S = zeros(size(X),'like',X);
    S(X > tol) = 1;
    S(X < -tol) = -1;
end