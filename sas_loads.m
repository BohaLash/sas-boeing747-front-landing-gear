
angle = linspace(0, pi/2, 100);
velocity = linspace(50, 150, 100);
[X, Y] = meshgrid(angle, velocity);
Q = force(X, Y);
surf(X,Y,Q)
xlabel("angle")
ylabel("velocity")
zlabel("force")

function P = weight() % N
    g = 9.8; % m/s^2
    m = 1432; % kg
    P = m * g; % N
end

function Mw = weight_moment(alpha) % Nm
    F = weight(); % N
    d0 = 0.5; % m
    Mw = F*d0*sin(alpha); % Nm
end

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
    D = 0.5*rho*s*Cd.*v.^2; % N
end

function L = lift(alpha, v) % N
    rho = 1.125;
    w = 0.5; % m
    h = 2; % m
    s = w*h; % m^2
    Cl = lift_coef(alpha);
    L = 0.5*rho*s*Cl.*v.^2; % N
end

function Md = drag_moment(alpha, v) % Nm
    D = drag(alpha, v); % N
    d0 = 1.5; % m
    Md = D*d0.*cos(alpha); % Nm
end

function Ml = lift_moment(alpha, v) % Nm
    L = lift(alpha, v); % N
    d0 = 1.5; % m
    Ml = -L*d0.*sin(alpha); % Nm
end

function l = dumper_length(alpha) % m
    b = 1; % m
    d = 0.5; % m
    m = 0.35; % m
    n1 = 0.45; % m
    n2 = 0.45; % m
    h = 0.1; % m
    gamma = acos((b^2 + d^2 - (n1 + n2)^2)/(2*b*d));
    
    e2 = b^2 + d^2 - 2*b*d*cos(gamma - alpha);
    e = sqrt(e2);
    eta = acos((e2 + n2^2 - n1^2)./(2*e*n2));
    zeta = asin(h/m);
    beta = eta - zeta;
    g2 = m^2 + h^2;
    g = sqrt(g2);
    f2 = g2 + e2 - 2*g*e.*cos(beta);
    f = sqrt(f2);

    l = real(f); % m
end

function Fd = dumper_force(alpha) % N
    K = 50000; % N/m
    f0 = 0.6; % m
    f = dumper_length(alpha); % m
    df = f0 - f; % m
    Fd = K*df; % N
end

function Mdp = dumper_moment(alpha) % N
    b = 1; % m
    d = 0.5; % m
    m = 0.35; % m
    n1 = 0.45; % m
    n2 = 0.45; % m
    h = 0.1; % m
    gamma = acos((b^2 + d^2 - (n1 + n2)^2)/(2*b*d));
    
    e2 = b^2 + d^2 - 2*b*d*cos(gamma - alpha);
    e = sqrt(e2);
    eta = acos((e2 + n2^2 - n1^2)./(2*e*n2));
    zeta = asin(h/m);
    beta = eta - zeta;
    g = sqrt(m^2 + h^2);
    f = dumper_length(alpha); % m
    omega = acos((e - g*cos(beta))./f);
    phi = acos((e2 + d^2 - b^2)./(2*e*d)) .* (1 - 2*(alpha < gamma));
    teta = pi - phi - omega;
    q = -sin(teta)*d; % m
    arm = real(q);
    F = dumper_force(alpha); % N
    Mdp = F.*arm; % Nm
end

function h = arm(alpha) % m
    a = 0.1; % m
    b = 1; % m
    gamma = pi/3; % rad
    h = a*b*sin(gamma + alpha)./sqrt(a^2 + b^2 - 2*a*b*cos(gamma + alpha)); % m
end

function f = force(alpha, v) % N
    M = weight_moment(alpha) + drag_moment(alpha, v) + lift_moment(alpha, v) + dumper_moment(alpha); % Nm
    % M = weight_moment(alpha) + dumper_moment(alpha); % Nm
    d = arm(alpha); % m
    f = M./d; % N
end
