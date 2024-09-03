function output = dilation_by_E(R,E)
    output = E * R * inv(E' * E) * E';
end