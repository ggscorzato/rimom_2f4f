function u = URI(mur0,Nf)

% Ciuchini et al Z. Phys. C 28 (1995) 239.
% simpler: Beciveric et. al EPJC 37 (2004) 315

Nc=3;
g0=4; 
b0=(11*Nc - 2*Nf)/3;
b1=(34/3)*Nc^2 - (10/3) * Nc*Nf - ((Nc^2-1)/Nc)*Nf;
LambMSr0=[0.60 0.62 ]; % hep-lat/0411025

logL = 2 * log(mur0./LambMSr0(Nf/2+1));
alpha= (4*pi)./(b0.*logL) .* (1 - (b1/b0^2) .* log(logL)./logL);
%alpha=alpha_s(mur0,3,LambMSr0,Nf,Nc)*4*pi;

J_RIMOM= 8 *log(2) - (17397 - 2070*Nf + 104* Nf^2)/(6*(33-2*Nf)^2);

u = alpha.^(-g0/(2*b0)) .* (1 + alpha/(4*pi) * J_RIMOM);