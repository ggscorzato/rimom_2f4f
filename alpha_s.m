function alpha_f = alpha_s(q,nl,Lambda,Nf,Nc)
% computes alpha_strong/(4*pi) at the scale q (which should be in the same unit as Lambda.
% nl = # of loop of precision % Nf, Nc: usual meaning

Z3=1.202056903;
Z4=1.082323234;
Z5=1.036927755;
C_f = (Nc^2-1)/(2*Nc);
%LambMSr0=[0.60 0.62];% hep-lat/0411025

b0  = 11/3*Nc - 2/3*Nf;
be1 = 34/3*Nc^2 - 38/3*Nf;
be2 = 2857/54*Nc^3 + C_f^2*Nf - 205/18*C_f*Nc*Nf -...
     1415/54*Nc^2*Nf + 11/9*C_f*Nf^2 + 79/54*Nc*Nf^2;

be3 = (150653/486 - 44/9*Z3)*Nc^4 +...
     (-39143/162 + 68/3*Z3)*Nc^3*Nf +...
     (7073/486 - 328/9*Z3)*C_f*Nc^2*Nf +...
     (-2102/27 + 176/9*Z3)*C_f^2*Nc*Nf +...
     23*C_f^3*Nf + (3965/162 + 56/9*Z3)*Nc^2*Nf^2 +...
     (338/27 - 176/9*Z3)*C_f^2*Nf^2 +...
     (4288/243 + 112/9*Z3)*C_f*Nc*Nf^2 + 53/243*Nc*Nf^3 +...
     154/243*C_f*Nf^3 +...
     (-10/27 + 88/9*Z3)*Nc^2*(Nc^2+36) +...
     (32/27 - 104/9*Z3)*Nc*(Nc^2+6)*Nf +...
     (-22/27 + 16/9*Z3)*(Nc^4 - 6*Nc^2 + 18)/Nc^2*Nf^2;

L2  = log(q.^2./Lambda.^2);
LL2 = log(L2);

b1 = be1/(b0*4*pi);
b2 = be2/(b0*16*pi^2);
b3 = be3/(b0*64*pi^3);

als0 = 4*pi./(b0.*L2);
als1 = als0 - als0.^2.*b1.*LL2;
als2 = als1 + als0.^3.*(b1^2*(LL2.^2 - LL2 -1) + b2);
als3 = als2 + als0.^4.*(b1^3*(-LL2.^3+5/2*LL2.^2+2*LL2-1/2)-3*b1*b2*LL2 + b3/2);

if( nl == 0 )
  alpha_f = als0/(4*pi);     %  Formula LO
elseif( nl == 1 ) 
  alpha_f = als1/(4*pi);     %  Formula NLO
elseif( nl == 2 ) 
  alpha_f = als2/(4*pi);     %  Formula N2LO
elseif( nl == 3 ) 
  alpha_f = als3/(4*pi);     %  Formula N3LO
end
