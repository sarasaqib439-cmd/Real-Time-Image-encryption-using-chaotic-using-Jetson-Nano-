function [Xout, Yout, Zout] = lorenz_rossler_system(delta, r, a, beta, b, c, X0, N)

dt = 0.01;
X = X0(:);

Xout = zeros(N,1);
Yout = zeros(N,1);
Zout = zeros(N,1);

for i = 1:N
    k1 = f(X);
    k2 = f(X + 0.5*dt*k1);
    k3 = f(X + 0.5*dt*k2);
    k4 = f(X + dt*k3);

    X = X + dt*(k1 + 2*k2 + 2*k3 + k4)/6;

    Xout(i) = X(1);
    Yout(i) = X(2);
    Zout(i) = X(3);
end

    function d = f(x)
        xx = x(1); yy = x(2); zz = x(3);

        dx = (delta-1)*yy - delta*xx - zz;
        dy = (r+1)*xx - (1-a)*yy - 20*xx*zz;
        dz = 5*xx*yy - beta*zz + b + xx*zz - c*xx;

        d = [dx; dy; dz];
    end
end
