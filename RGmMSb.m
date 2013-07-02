function m_o = RGmMSb(mu_ini,mu_fin, Lambda, m_i, Nf)
% RG evolution for the mass, MSbar

m_o = (cm(mu_fin,Lambda,Nf)./cm(mu_ini,Lambda,Nf)) .*m_i;

function c= cm(q,Lambda,Nf)
nl=3;
Nc=3;
al = alpha_s(q, nl, Lambda, Nf,Nc)*4*pi;

if(Nf==0)
c=al^(4/11) * (1 + 0.218783 * al + 0.153208 * al^2 + 0.130873 * al^3);
elseif(Nf==2)
c=al^(12/29) * (1 + 0.256553 * al + 0.141946 * al^2 + 0.0880197 * al^3);
elseif(Nf==3)
c=al^(4/9) * (1 + 0.284907 * al + 0.138955 * al^2 + 0.0629448 * al^3);
elseif(Nf==4)
c=al^(12/25) * (1 + 0.322809* al + 0.140756 * al^2 + 0.0351717 * al^3);
end