function m_o = RGmRI(mu_ini,mu_fin, Lambda, m_i, Nf)
% RG evolution for the mass, RI

m_o = (cm(mu_fin,Lambda,Nf)./cm(mu_ini,Lambda,Nf)) .*m_i;


function c= cm(q,Lambda,Nf)
nl=3;
Nc=3;
al = alpha_s(q, nl, Lambda, Nf,Nc)*4*pi;

if(Nf==0)
c=al^(4/11) * (1 + 0.643196 * al + 1.44071 * al^2 + 4.46061 * al^3);
elseif(Nf==2)
c=al^(12/29) * (1 + 0.680966 * al + 1.32023 * al^2 + 3.55806 * al^3);
elseif(Nf==3)
c=al^(4/9) * (1 + 0.70932 * al + 1.26666 * al^2 + 3.13088 * al^3);
elseif(Nf==4)
c=al^(12/25) * (1 + 0.747222 * al + 1.22192 * al^2 + 2.72228 * al^3);
end