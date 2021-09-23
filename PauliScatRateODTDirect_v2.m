function [outODT] = PauliScatRateODTDirect_v2(T, kappa, N, truncationParam)
% Normalized scattering rate due to Pauli blocking in an ODT potential, based on the local density approximation
% We integrate the box potential Pauli scattering rate which is a function of density, where the weight function is the ODT density
% TTF = T/TF
% kappa is the dimsionless momentum
% Function cannot accept vectors!!
% Use the truncation of the potential, where truncationParam determines how many kB*T to truncate from the potential
% v2: use truncation by setting potential walls to infinity: UODT(UODT > (U0-deltaU)) = Inf;

kB = 1.3806504e-23; %J/K
hbar = 1.054571628e-34;
m = 6/6.022e26;
g = 2*pi*5.8724e6;
lambdaAtom = 670.977338e-9;
c = 2.99792458e8; % speed of light [m/s]

Vcontrol = 0.5; % ODT control voltage

Trans = 0.9936 * 0.999* 0.9926 * 0.995 * 0.9975; %Transmission of two lenses (AC508-500-C-ML and AC508-250-C-ML), one laser line mirror, and two surfaces of the vacuum window (from Aviv's thesis, Appendix C, p. 125).
P = Trans.*(1.7893.*Vcontrol + 0.0014); % taking into account optical transmissions [W]
w0 = 7.9786e-6;
lambdaIR = 1064e-9;
omega0 = 2*pi*c/lambdaAtom;
omega = 2*pi*c/lambdaIR;
zR = pi*w0*w0/lambdaIR * 2684/2019; %the extra factor is needed in order to explain the low measured axial frequency (2020-02-27)
U0 = 3*pi*c^2/(2*omega0^3)*( g/(omega0-omega) + g/(omega0+omega))*2.*P./(pi*w0*w0);

omega_r = sqrt(4.*U0./(m*w0*w0));
omega_z = sqrt(2.*U0./(m*zR^2));

Ef = hbar .* (omega_z .* omega_r .* omega_r .* 6 .* N).^(1/3); %Fermi energy [J]
Tf = Ef./kB; % Fermi temperature [K]
TTF = T./Tf;

deltaU = truncationParam.*kB.*T;
% deltaU = 0.02*U0; % override for 98% of potential depth
muODT = FugacityAnharmonicInternal(T, N, deltaU);
zODT = exp(muODT./(kB.*T));
% zHarmonic = fugacity(TTF);
% muHarmonic = qHarmonic.*kB*T;


%% Integration limits
% Integration ranges: a (0, Inf); y (-Inf, Inf), r (0, IntSizeR*w0), z (-IntSizeZ*zR,IntSizeZ*zR)
[rmax, zmax] = integrationLimits(U0, w0, zR, deltaU);

%% Harmonic case - to test the code with respect to the known result:
% intergrandHarmonic  = @(a, y, r, z) r ./ (1 + zHarmonic.^(-1).*exp( a + y.^2 + 1./(kB*T).*UHarmonic(U0, w0, r, z, zR)) )  .* 1./(1 + zHarmonic.^(-1).*exp( a + (y + kappa./sqrt(TTF)).^2 + 1./(kB*T).*UHarmonic(U0, w0, r, z, zR) ) );
% Seems to work fine - agrees with Thywissen's result
% limit = Inf;
% intHarmonic1 = @(z) integral3(@(a, y, r) intergrandHarmonic(a, y, r, z), 0, limit, -limit, limit,0, IntSizeR*w0,'AbsTol', 1e-3,'RelTol',1e-10, 'method', 'iterated'); %
% intHarmonic2 = integral(intHarmonic1, -IntSizeZ*zR,IntSizeZ*zR , 'ArrayValued', true, 'AbsTol', 1e-3,'RelTol',1e-10);
% disp(['limit = ' num2str(limit) ', int = ' num2str(1 - (2*m*kB*T)^(3/2)*(2*pi)^2./( 2 .* (2*pi*hbar)^3.*N ) .* intHarmonic2)])
% outHarmonic = 1 - (2*m*kB*T)^(3/2)*(2*pi)^2./( 2 .* (2*pi*hbar)^3.*N ) .* intHarmonic2;

%% ODT case
intergrandODT = @(a, y, r, z) r ./ (1 + zODT.^(-1).*exp( a + y.^2 + 1./(kB*T).*UODT(U0, w0, r, z, zR)) )  .* 1./(1 + zODT.^(-1).*exp( a + (y + kappa./sqrt(TTF)).^2 + 1./(kB*T).*UODT(U0, w0, r, z, zR) ) );

intODT1 = @(z) integral3(@(a, y, r) intergrandODT(a, y, r, z), 0, Inf, -Inf, Inf, 0, rmax,'AbsTol', 1e-3,'RelTol',1e-10, 'method', 'iterated');
intODT2 = integral(intODT1, -zmax, zmax, 'ArrayValued', true, 'AbsTol', 1e-3,'RelTol',1e-10);

outODT = 1 - (2*m*kB*T)^(3/2)*(2*pi)^2./( 2 .* (2*pi*hbar)^3.*N ) .* intODT2;
disp(['Gamma ODT = ' num2str(outODT)])


%% test integrals - should give 1 if normalization is correct (integration limits should be the same as in fugacity calculation - 95% or 98%)
% testIntergrandHarmonic = @(r,z) 2*pi.*r .* 1 .* nFHarmonic(r,z,T,U0);
% testIntergrandODT = @(r,z) 2*pi.*r .* 1 .* nFDOT(r,z,T,U0);
% testIntHarmonic = 1 - 1./N .* integral2(testIntergrandHarmonic, 0,IntSizeR*w0,-IntSizeZ*zR,IntSizeZ*zR,'AbsTol', 1e-9,'RelTol',1e-9, 'method', 'iterated');
% testIntODT = 1 - 1./N .* integral2(testIntergrandODT, 0,IntSizeR*w0,-IntSizeZ*zR,IntSizeZ*zR,'AbsTol', 1e-9,'RelTol',1e-9, 'method', 'iterated');
% disp(['Harmonic error: ' num2str(testIntHarmonic)])
% disp(['ODT error: ' num2str(testIntODT)])

%%
% [R,Z] = meshgrid(linspace(0,0.5*w0,20),linspace(-zR,zR,20));
% 
% figure;surf(R,Z,intergrandHarmonic(R,Z))

end

function UODT = UODT(U0, w0, r, z, zR) %deltaU
% ODT potential
% deltaU is the truncation parameter, such that UODT(r,z)<=U0-deltaU
UODT = U0 - U0./( 1 + (z./zR).^2 ) .* exp(-2.*r.^2./ (w0^2*( 1 + (z./zR).^2 )) );
% UODT(UODT > (U0-deltaU)) = Inf; %Truncate potential
end

function UHarmonic = UHarmonic(U0, w0, r, z, zR)
% Harmonic approximation of the ODT potential:
UHarmonic = U0 - U0.*( 1 - (z./zR).^2 - 2.*(r./w0).^2 );
end

function [muOut, mu0] = FugacityAnharmonicInternal(T, N, deltaU)
% T is the temperautre
% N is is the required / requested number of atoms, for which we ned to calculate the checmical potential mu
% Important: to agree with regular convenstions (e.g. Varenna notes), potentials are defined such that the energy at z=0, r=0 is 0 (!!!)
% Output parameters: muOut is the calculated chemical potential; mu0is the harmonic chemical potential
kB = 1.3806504e-23; %J/K
hbar = 1.054571628e-34;
% a0 = 0.52917720859e-10; %Bohr radius [m]
m = 6/6.022e26;
g = 2*pi*5.8724e6;
lambdaAtom = 670.977338e-9;
c = 2.99792458e8; % speed of light [m/s]
% mu0 = 1.25663706212e-6; %Vacuum permeability [H/m]
% muB = 9.27400915e-24; %Bohr magneton, J/T
% lambdaDB = sqrt(2*pi*hbar^2./ (m.*kB.*T) );

Vcontrol = 0.5; % ODT control voltage

Trans = 0.9936 * 0.999* 0.9926 * 0.995 * 0.9975; %Transmission of two lenses (AC508-500-C-ML and AC508-250-C-ML), one laser line mirror, and two surfaces of the vacuum window (from Aviv's thesis, Appendix C, p. 125).
P = Trans.*(1.7893.*Vcontrol + 0.0014); % taking into account optical transmissions [W]
w0 = 7.9786e-6;
lambdaIR = 1064e-9;
omega0 = 2*pi*c/lambdaAtom;
omega = 2*pi*c/lambdaIR;
zR = pi*w0*w0/lambdaIR * 2684/2019; %the extra factor is needed in order to explain the low measured axial frequency (2020-02-27)
U0 = 3*pi*c^2/(2*omega0^3)*( g/(omega0-omega) + g/(omega0+omega))*2.*P./(pi*w0*w0);

omega_r = sqrt(4.*U0./(m*w0*w0));
omega_z = sqrt(2.*U0./(m*zR^2));

Ef = hbar .* (omega_z .* omega_r .* omega_r .* 6 .* N).^(1/3); %Fermi energy [J]
Tf = Ef./kB; % Fermi temperature [K]

[~, q] = fugacity(T./Tf); %Harmonic value
mu0 = q.*kB*T; % Initial guess of the chemical potential, based on the harmonic approximation, when using potential definition of U(r=0,z=0)=0

%integration limits on r and z, up to 98% of the potential - to avoid integrating outside of the trap:
% IntSizeR = 1.39857;
% IntSizeZ = 7;
%integration limits on r and z, up to 95%
% IntSizeR = 1.22387;
% IntSizeZ = 4.3589;

[rmax, zmax] = integrationLimits(U0, w0, zR, deltaU); %Calculate integration limits (box shape), according to U(rmax,z=0)=U0-deltaU, and U(r=0,zmax)=U0-deltaU

nFODT = @(r,z,T,U0,mu) -2*pi*r./(2*pi*hbar^2./ (m.*kB.*T) ).^(3/2) .* polylog3half(- exp( 1./(kB*T) .*(mu - UODT(U0, w0, r, z, zR)) ));
% AtomNo = integral2( @(r,z)nFODT(r,z,T,U0), 0,IntSizeR*w0,-IntSizeZ*zR,IntSizeZ*zR,'AbsTol', 1e-20,'RelTol',1e-20); %,'AbsTol', 1e-20,'RelTol',1e-20
AtomNoIntegral = @(mu) integral2( @(r,z)nFODT(r,z,T,U0, mu), 0, rmax, -zmax, zmax, 'AbsTol', 0.1,'RelTol',1e-8, 'method', 'iterated'); %,'AbsTol', 1e-20,'RelTol',1e-20

AtomNumDiffAbs = @(mu) abs(N - AtomNoIntegral(mu));
options = optimset('MaxFunEvals', 1e6, 'MaxIter', 1e6);
muOut = fminsearch( @(mu)AtomNumDiffAbs(mu), mu0, options);


end

function [rmaxOut, zmaxOut] = integrationLimits(U0, w0, zR, deltaU)
% Find rmax and zmax, such that U(rmax,z=0)=U0-deltaU, and U(r=0,zmax)=U0-deltaU (this leads to box integration limits)
syms rmax
syms zmax

rmaxOut = abs(double( vpasolve( UODT(U0, w0, rmax, 0, zR)==U0-deltaU, rmax, w0 ) ));
zmaxOut = abs(double( vpasolve( UODT(U0, w0, 0, zmax, zR)==U0-deltaU, zmax, zR ) ));
if length(rmaxOut)>1
    rmaxOut = rmaxOut(1);
end
if length(zmaxOut)>1
    zmaxOut = zmaxOut(1);
end
end