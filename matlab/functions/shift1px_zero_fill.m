function Y = shift1px_zero_fill(X, dRow, dCol)
    Y = false(size(X));

    rSrc = (1:size(X,1)) - dRow;
    cSrc = (1:size(X,2)) - dCol;

    rOK = rSrc>=1 & rSrc<=size(X,1);
    cOK = cSrc>=1 & cSrc<=size(X,2);

    Y(rOK, cOK) = X(rSrc(rOK), cSrc(cOK));
end