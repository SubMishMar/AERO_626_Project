function Z = measurement_model(X)
x = X(1);
y = X(2);
z = X(3);

Z = x^2 + y^2;
end