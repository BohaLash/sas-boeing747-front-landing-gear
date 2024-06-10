% rho = 1.125;
% w = 0.5; % m
% h = 2; % m
% s = w*h; % m^2
% v = 55; % m/s
% Cd = 0.15;
% D = 0.5*rho*s*Cd*v^2; % N
% d0 = 1; % m
% M = D*d0; % Nm
% 
% [D, M]

% a = 0.1;
% gamma = 1.57;
% d = a * cos(gamma - 0.52)

cfd_data = readtable('sas_aero.csv');
cfd_lookup_alpha = (90 - cfd_data.a) / 360 * 2 * pi;
cfd_lookup_cd = cfd_data.Cd;
cfd_lookup_cl = cfd_data.Cl;

v = 60; %m/s

force(0, v)
force(1.57, v)

x = linspace(0, 1.57, 100);
y = 1:100;
for i = 1:100
    % y(i) = interp1(cfd_lookup_alpha, cfd_lookup_cd, x(i));
    % y(i) = interp1(cfd_lookup_alpha, cfd_lookup_cl, x(i));
    % y(i) = drag(x(i), v);
    % y(i) = arm(x(i));
    y(i) = force(x(i), v);
end

plot(x, y)

function Cd = drag_coef(alpha) % N
    cfd_data = readtable('sas_aero.csv');
    cfd_lookup_alpha = (90 - cfd_data.a) / 360 * 2 * pi;
    cfd_lookup_cd = cfd_data.Cd;
    Cd = interp1(cfd_lookup_alpha, cfd_lookup_cd, alpha);
end

function Cl = lift_coef(alpha) % N
    cfd_data = readtable('sas_aero.csv');
    cfd_lookup_alpha = (90 - cfd_data.a) / 360 * 2 * pi;
    cfd_lookup_cl = cfd_data.Cl;
    Cl = interp1(cfd_lookup_alpha, cfd_lookup_cl, alpha);
end

function D = drag(alpha, v) % N
    rho = 1.125;
    w = 0.5; % m
    h = 2; % m
    s = w*h; % m^2
    Cd = drag_coef(alpha);
    D = 0.5*rho*s*Cd*v^2; % N
end

function L = lift(alpha, v) % N
    rho = 1.125;
    w = 0.5; % m
    h = 2; % m
    s = w*h; % m^2
    Cl = lift_coef(alpha);
    L = 0.5*rho*s*Cl*v^2; % N
end

function Md = drag_moment(alpha, v) % Nm
    D = drag(alpha, v); % N
    d0 = 1.5; % m
    Md = D*d0*cos(alpha); % Nm
end

function Ml = lift_moment(alpha, v) % Nm
    L = lift(alpha, v); % N
    d0 = 1.5; % m
    Ml = -L*d0*sin(alpha); % Nm
end

function h = arm(alpha) % m
    a = 0.1; % m
    b = 1; % m
    gamma = pi/3; % rad
    h = a*b*sin(gamma + alpha)./sqrt(a^2 + b^2 - 2*a*b*cos(gamma + alpha)); % m
end

function f = force(alpha, v) % N
    M = lift_moment(alpha, v) + drag_moment(alpha, v); % Nm
    d = arm(alpha); % m
    f = M/d; % N
end
