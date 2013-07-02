function [go,dgo,Zo,dZo] = computeZ2(prop,L,gam_frmt,Zpsi,dZpsi)
%gam_frmt = 'milc', 'gwc ', 
addpath /Users/luigi/Work/lib/Matlab

D=length(L);
V=prod(L);
sV=prod(L(1:D-1));
NS=2^(D/2);
T=L(D);
NC=3;
CS=NC*NS;

%%%%%%%%%%%%%% gamma matrices
[ga gi gv] = gammat(gam_frmt);

%%%%%%%%%%%%%% compute Momenta
[P,Pb,Ph,P2,Pb2,Ph2]= computeMom(L);

%%%%%%%%%%%%%%% initialize

nconf=size(prop,6);
maxmom=size(prop,5);

prop=permute(prop,[1 3 2 4 5 6]);  % c_sink s_sink c_src s_src 
prop=reshape(prop,CS,CS,maxmom,nconf);
G=zeros(CS,CS,11);
avProp=zeros(CS,CS,maxmom);
avProp=mean(prop,4);

%%%%%%%%%%%%%% compute Zpsi if not aready available

if(isempty(Zpsi))
  nbootsamp=200;
  [Zpsi,dZpsi] = computeZpsi(prop,L,gam_frmt,nbootsamp);
end

%%%%%%%%%%%%% loop over momenta
for ip=1:maxmom
  %% S(p) and S^-1(p)
  invProp(:,:) =  inv(avProp(:,:,ip));

  %% G
  for n=1:nconf
    G(:,:,1) =  (prop(:,:,ip,n) * ga(:,:,5) * ctranspose(prop(:,:,ip,n)) * ga(:,:,5));
    G(:,:,2) =  (prop(:,:,ip,n) * ctranspose(prop(:,:,ip,n)) * ga(:,:,5));
    for o=1:4
      G(:,:,2+o) = (prop(:,:,ip,n) * ga(:,:,o) * ga(:,:,5) * ctranspose(prop(:,:,ip,n)) * ga(:,:,5));
      G(:,:,6+o) = (prop(:,:,ip,n) * ga(:,:,o) * ctranspose(prop(:,:,ip,n)) * ga(:,:,5));
    end
    G(:,:,11) = (prop(:,:,ip,n) * ctranspose(prop(:,:,ip,n)));

    %% Lambda
    for o=1:11
      Lambda(:,:,o)=invProp(:,:) * G(:,:,o) * invProp(:,:); 
    end

    %% Gamma
    Gamma(ip,1,n) = real(trace(Lambda(:,:,1))) /CS; % S
    Gamma(ip,2,n) = real(trace(Lambda(:,:,2) * ga(:,:,5))) /CS; % P
    for o=1:4
      Gamma(ip,2+o,n) = real(trace(Lambda(:,:,2+o) * ga(:,:,o))) /(CS*4); % V
      Gamma(ip,6+o,n) = real(trace(Lambda(:,:,6+o) * ga(:,:,5) * ga(:,:,o))) /(CS*4); % A
    end
    Gamma(ip,11,n) = real(trace(Lambda(:,:,11))) /CS;
    Gamma(ip,12,n) = real(trace(G(:,:,11))) /CS;
    
  end %conf

  %% Gam
  Gam(ip,1)=mean(squeeze(Gamma(ip,1,:)));
  Gam(ip,2)=mean(squeeze(Gamma(ip,2,:)));
  Gam(ip,3)=mean(squeeze(sum(Gamma(ip,[3: 6],:),2)));
  Gam(ip,4)=mean(squeeze(sum(Gamma(ip,[7:10],:),2)));

  dGam(ip,1)=std(squeeze(Gamma(ip,1,:)))./sqrt(nconf);
  dGam(ip,2)=std(squeeze(Gamma(ip,2,:)))./sqrt(nconf);
  dGam(ip,3)=std(squeeze(sum(Gamma(ip,[3: 6],:),2)))./sqrt(nconf);
  dGam(ip,4)=std(squeeze(sum(Gamma(ip,[7:10],:),2)))./sqrt(nconf);

end %ip

for o=1:4
  [go(:,o) dgo(:,o)]=mean_over_equalP2(Gam(:,o),dGam(:,o),P2);
end

  %%% Z_o
  nf=2;
for ip=1:maxmom
  for o=1:4 
    Zo(ip,o)= Zpsi(ip)^(nf/2) / go(ip,o);
    dZo(ip,o)= Zo(ip,o)*sqrt((dgo(ip,o)/go(ip,o))^2 + (dZpsi(ip)/Zpsi(ip))^2);
  end
end % ip
