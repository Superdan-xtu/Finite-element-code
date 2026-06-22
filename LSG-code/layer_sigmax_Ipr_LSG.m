function [s11x, s22x, s12x] = layer_sigmax_Ipr_LSG(X, La, mu)

x = X(:,1);
y = X(:,2);

S = sin(1);
E = exp(1);

sx = sin(x);  cx = cos(x);  ex = exp(x);
sy = sin(y);  cy = cos(y);  ey = exp(y);

% u1 = F(x) G(y)
F  = sx .* (sx - S);
G  = sy .* (sy - S);

Fp  = cx .* (2*sx - S);
G_p = cy .* (2*sy - S);

Fpp = 2*cx.^2 - 2*sx.^2 + S*sx;

% u2 = H(x) K(y)
H  = sx .* (ex - E);
K  = sy .* (ey - E);

Hp  = cx .* (ex - E) + sx .* ex;
Kp  = cy .* (ey - E) + sy .* ey;

Hpp = E*sx + 2*ex.*cx;

% second derivatives needed for x-derivatives of sigma
u1xx = Fpp .* G;
u1xy = Fp  .* G_p;
u2xx = Hpp .* K;
u2xy = Hp  .* Kp;

% sigma11 = (2mu+lambda)u1x + lambda u2y
% sigma22 = lambda u1x + (2mu+lambda)u2y
% sigma12 = mu(u1y+u2x)
s11x = (2*mu + La)*u1xx + La*u2xy;
s22x = La*u1xx + (2*mu + La)*u2xy;
s12x = mu*(u1xy + u2xx);

end
