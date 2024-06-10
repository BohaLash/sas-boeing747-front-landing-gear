clc; clear; close all;

a = 0.1;
b = 1;
d = 0.5;
m = 0.35;
n = 0.45;
n1 = 0.45;
n2 = 0.45;
h = 0.1;
gamma = acos((b^2 + d^2 - (n1 + n2)^2)/(2*b*d));

gamma_deg = gamma * 360 / 2 / pi

alpha = linspace(0, pi/2, 100000);

e2 = b^2 + d^2 - 2*b*d*cos(gamma - alpha);
e = sqrt(e2);
eta = acos((e2 + n2^2 - n1^2)./(2*e*n2));
zeta = asin(h/m);
beta = eta - zeta;
g2 = m^2 + h^2;
g = sqrt(g2);
f2 = g2 + e2 - 2*g*e.*cos(beta);

f = sqrt(f2); % m

% omega = acos((e2 + f2 - g2)./(2*e.*f));
omega = acos((e - g*cos(beta))./f);
% phi = acos((e2 + d^2 - b^2)./(2*e*d));
phi = acos((e2 + d^2 - b^2)./(2*e*d)) .* (1 - 2*(alpha < gamma));
teta = pi - phi - omega;

q = -sin(teta)*d; % m
f0 = 0.6; % m
K = 1; % N/m
F = K*(f0 - f); % N
M = F.*q; % Nm

x = alpha * 360 / 2 / pi;
y = M;
plot(x, y)
