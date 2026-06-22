function [s11y, s22y, s12y] = layer_sigmay_Ipr_LSG(X, La, mu)

x = X(:,1);
y = X(:,2);

S = sin(1);
E = exp(1);

sx = sin(x);  cx = cos(x);  ex = exp(x);
sy = sin(y);  cy = cos(y);  ey = exp(y);

% u1 = F(x) G(y)
F  = sx .* (sx - S);
G  = sy .* (sy - S);

Fp = cx .* (2*sx - S);
Gp = cy .* (2*sy - S);

Gpp = 2*cy.^2 - 2*sy.^2 + S*sy;

% u2 = H(x) K(y)
H  = sx .* (ex - E);
K  = sy .* (ey - E);

Hp = cx .* (ex - E) + sx .* ex;
Kp = cy .* (ey - E) + sy .* ey;

Kpp = E*sy + 2*ey.*cy;

% second derivatives needed for y-derivatives of sigma
u1xy = Fp .* Gp;
u1yy = F  .* Gpp;
u2xy = Hp .* Kp;
u2yy = H  .* Kpp;

% sigma11 = (2mu+lambda)u1x + lambda u2y
% sigma22 = lambda u1x + (2mu+lambda)u2y
% sigma12 = mu(u1y+u2x)
s11y = (2*mu + La)*u1xy + La*u2yy;
s22y = La*u1xy + (2*mu + La)*u2yy;
s12y = mu*(u1yy + u2xy);

end
