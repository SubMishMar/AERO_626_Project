function Z = measurement_model(X)
x = X(1);
y = X(2);
z = X(3);

Z = [x; y];
end