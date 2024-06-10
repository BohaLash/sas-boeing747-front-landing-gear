0.1 + 1 - sqrt(0.1^2 + 1^2 - 2*0.1*1*cos(pi/3))

[gear_angle(0) gear_angle(0.146)]

x = linspace(0, 0.146, 1000);
y = gear_angle(x);

plot(x, y)

function ang = gear_angle(x)
    a = 0.1;
    b = 1;
    gamma = pi/3;
    c = sqrt(a*a + b*b - 2*a*b*cos(gamma));
    ang = acos((a.^2 + b.^2 - (c + x).^2)/(2 * a * b)) - gamma;
end
