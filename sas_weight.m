% g = 9.8; % m/s^2
% m = 1432; % kg
% F = m * g; % N
% d0 = 1; % m
% M = F*d0; % Nm

% a = 0.1;
% gamma = 1.57;
% d = a * cos(gamma - 0.52)

force(0)
force(1.57)

x = linspace(0, 1.57, 100);
y = 1:100;
for i = 1:100
    y(i) = force(x(i));
end

plot(x, y)

function P = weight() % N
    g = 9.8; % m/s^2
    m = 1432; % kg
    P = m * g; % N
end

function M = moment(alpha) % Nm
    F = weight(); % N
    d0 = 0.5; % m
    M = F*d0*sin(alpha); % Nm
end

function h = arm(alpha) % m
    a = 0.1; % m
    b = 1; % m
    gamma = pi/3; % rad
    h = a*b*sin(gamma + alpha)./sqrt(a^2 + b^2 - 2*a*b*cos(gamma + alpha)); % m
end

function f = force(alpha) % N
    M = moment(alpha); % Nm
    d = arm(alpha); % m
    f = M/d; % N
end
