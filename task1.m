%%%task 1
g = sym('g','real'); %gamma
c = [0; 2*g; 1; 1];
bHat = [(-4*g^2+6*g-1)/(4*g); (-2*g+1)/(4*g); g; 0];
b = [(6*g-1)/(12*g); -1/(12*g*(2*g-1)); (-6*g^2+6*g-1)/(3*(2*g-1)); g];
A = [0,0,0,0;...
    g,g,0,0;...
    (-4*g^2+6*g-1)/(4*g), (-2*g+1)/(4*g), g, 0;...
    (6*g-1)/(12*g), -1/(12*g*(2*g-1)), (-6*g^2+6*g-1)/(3*(2*g-1)), g];

fprintf('\n\n');
fprintf('test order 1:\n\n');
fprintf('(sum(b_i) = 1):\n');
fprintf(' - b    => %d\n',simplify(sum(b)));
fprintf(' - bHat => %d\n',simplify( sum(bHat)));

fprintf('\n\n');
fprintf('test order 2:\n\n');
fprintf('(sum(b_i*c_i) = 1/2):\n')
fprintf(' - b    => %s\n',simplify(sum(b.*c)));
fprintf(' - bHat => %s\n',simplify(sum(bHat.*c)));

fprintf('\n\n')
fprintf('test order 3:\n\n')
fprintf('(sum(b_i*c_i^2) = 1/3):\n')

fprintf(' - b    => %s\n',simplify(sum(b.*c.^2)));
fprintf(' - bHat => ');
disp(simplify(sum(bHat.*c.^2)));
fprintf('Error estimating method is of order 2 for "gamma" != 1/2 +- 1/(2*sqrt(3)).\n')

fprintf('\n\n');
fprintf('(sum(b_i*a_ij*c_j) = 1/6):\n');
fprintf(' - b    => %s\n',simplify(sum(b'*A*c)));

fprintf('\n\n')
fprintf('test order 4:\n\n')
fprintf('(sum(b_i*c_i^3) = 1/4):\n')
fprintf(' - b    => %s\n',simplify(sum(b.*c.^3)));
fprintf('Advancing method is of order 3 for "gamma" != 1/4.\n')