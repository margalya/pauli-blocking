function [out] = PauliScatRate(TTF, kappa)
% Normalized scattering rate due to Pauli blocking in an harmonic potential
% Integral result taken from: B. Shuve and J. H. Thywissen, J. Phys. B At. Mol. Opt. Phys. 43, 015301 (2010).
% TTF is T/TF
% kappa is the dimensionless momentum
% Using the explicit integral equation for the fugacity instead of simply the Sommefeld approximation
% fun1 = @(a, y, TTF, kappa) TTF.^3./(1 + exp(a + y.^2 - (1 - pi.^2./3.*TTF.^2)./TTF)).*a.^(3/2)./( 1 + exp(-(1 - pi.^2./3.*TTF.^2)./TTF + a + (y + kappa./sqrt(TTF)).^2));

% Beyond Sommerfeld - exact equation for z:
z = fugacityInternal(TTF);
fun1 = @(a, y, TTF, kappa, z) TTF.^3 .* a.^(3/2) ./ (1 + z.^(-1).*exp( a + y.^2) )  .* 1./(1 + z.^(-1).*exp( a + (y + kappa./sqrt(TTF)).^2 ) );

% q = integral3(fun,xmin,xmax,ymin,ymax,zmin,zmax)
int = zeros(size(TTF));
if any(size(kappa)>1) %kappa is an array
    for i = 1 : length(TTF)
        int(i) = integral2(@(a,y)fun1(a, y, TTF(i), kappa(i), z(i)), 0,inf ,-inf,inf,'AbsTol', 1e-20,'RelTol',1e-20);
    end
else
    for i = 1 : length(TTF)
        int(i) = integral2(@(a,y)fun1(a, y, TTF(i), kappa, z(i)), 0,inf ,-inf,inf,'AbsTol', 1e-20,'RelTol',1e-20);
    end
end
out = 1 - 8./pi*int;

end

function [z, q] = fugacityInternal(TTF)
% Calculate the fugacity z=exp(m*beta) and the logarithmic fugacity q=mu*beta of a gas as a function of T/TF = TTF
syms z
fugacity = zeros(size(TTF));
for i = 1:length(TTF)
fugacity(i) = double(vpasolve(polylog(3,-z)+1./6./(TTF(i)).^3));
end
clear z
z = fugacity;
q = log(fugacity);

end
