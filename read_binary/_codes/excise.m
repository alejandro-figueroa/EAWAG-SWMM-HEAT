% removes NaN values

function X = excise(X)
X(any(isnan(X)'),:) = [];