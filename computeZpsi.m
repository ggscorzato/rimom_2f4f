function [Zpsi,dZpsi] = computeZpsi(prop,L,gam_frmt,nbootsamp)
%gam_frmt = 'milc', 'gwc ', 
addpath /Users/luigi/Work/lib/Matlab

D=length(L);
V=prod(L);
sV=prod(L(1:D-1));
NS=2^(D/2);
T=L(D);
NC=3;
CS=NC*NS;
ISOVERLAP=0;

%%%%%%%%%%%%%% gamma matrices
[ga gi gv] = gammat(gam_frmt);

%%%%%%%%%%%%%% compute Momenta
[P,Pb,Ph,P2,Pb2,Ph2]= computeMom(L);
if(ISOVERLAP)
  rho=1;
  omega= sqrt(Pb2 + (Ph2/2 - rho).^2)./rho;
else
  omega=ones(maxmom);
end

%%%%%%%%%%%%%%% compute Z

nconf=size(prop,6);
maxmom=size(prop,5);
prop=permute(prop,[1 3 2 4 5 6]); % c_sink s_sink c_src s_src
prop=reshape(prop,CS,CS,maxmom,nconf);
avProp=zeros(CS,CS,maxmom);
avProp=mean(prop,4);
BavProp=zeros(CS,CS);
Zpsi=zeros(maxmom,1);
dZpsi=nan*ones(maxmom,1);

for ip=2:maxmom
  %% S(p) and S^-1(p)
  invProp(:,:) =  inv(avProp(:,:,ip));

  %%% Zpsi
  uu=zeros(CS,CS);
  for mu=1:4
    uu=uu+(Pb(mu,ip) * (ga(:,:,mu) * invProp(:,:)));
  end
  Zpsi(ip) = real(-i/CS * (trace(uu) /Pb2(ip))*omega(ip));
end


if(nbootsamp>0)
for ip=2:maxmom
  for ib=1:nbootsamp
    ii=ceil(nconf*rand(nconf,1));
    BavProp=mean(squeeze(prop(:,:,ip,ii)),3);

    %% S(p) and S^-1(p)
    invProp(:,:) =  inv(BavProp(:,:));

    %%% Zpsi
    uu=zeros(CS,CS);
    for mu=1:4
      uu=uu+(Pb(mu,ip) * (ga(:,:,mu) * invProp(:,:)));
    end
    Zb(ib) = real(-i/CS * (trace(uu) /Pb2(ip))*omega(ip));

  end %ib
  dZpsi(ip) = std(Zb);
end % ip
end % if nbootsamp > 0


[Zpsi dZpsi]=mean_over_equalP2(Zpsi,dZpsi,P2);

