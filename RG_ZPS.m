function Z_fin = RG_ZPS(mu_ini,mu_fin,nl, Lambda, Z_ini, Nf)
Nc=3;
alm = alpha_s(mu_fin, nl, Lambda, Nf,Nc);
al0 = alpha_s(mu_ini, nl, Lambda, Nf,Nc);

if(nl == 0)
  if( Nf == 0 )
    cmu = (alm)^(-4./11.);       % Nf = 0
    cm0 = (al0)^(-4./11.);
  elseif( Nf == 2 )
    cmu = (alm)^(-12./29.);      % Nf = 2
    cm0 = (al0)^(-12./29.);
  end
elseif(nl == 1)
  if( Nf == 0 ) 
    cmu = (alm)^(-4./11.) * (1. - 8.0826 * alm);
    cm0 = (al0)^(-4./11.) * (1. - 8.0826 * al0);
  elseif( Nf == 2 ) 
    cmu = (alm)^(-12./29.) * (1. - 8.55727 * alm);
    cm0 = (al0)^(-12./29.) * (1. - 8.55727 * al0);
  end
elseif(nl == 2) 
  if( Nf == 0 ) 
    cmu = (alm)^(-4./11.) * (1. - 8.0826 * alm - 162.179 * alm^2);
    cm0 = (al0)^(-4./11.) * (1. - 8.0826 * al0 - 162.179 * al0^2);
  elseif( Nf == 2 ) 
    cmu = (alm)^(-12./29.) * (1. - 8.55727 * alm - 135.256 * alm^2);
    cm0 = (al0)^(-12./29.) * (1. - 8.55727 * al0 - 135.256 * al0^2);
  end
elseif(nl == 3) 
  if( Nf == 0 ) 
    cmu = (alm)^(-4./11.) * (1. - 8.0826 * alm - 162.179 * alm^2 - 5701.94 * alm^3);
    cm0 = (al0)^(-4./11.) * (1. - 8.0826 * al0 - 162.179 * al0^2 - 5701.94 * al0^3);
  elseif( Nf == 2 ) 
    cmu = (alm)^(-12./29.) * (1. - 8.55727 * alm - 135.256 * alm^2 -4119.14 * alm^3);
    cm0 = (al0)^(-12./29.) * (1. - 8.55727 * al0 - 135.256 * al0^2 -4119.14 * al0^3);
  end
end

Z_fin = cm0/cmu*Z_ini;


