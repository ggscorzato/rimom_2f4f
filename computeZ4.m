function [Zpsi Z4F GAM_a, GAM_d] = computeZ4(prop,L,gam_frmt)
%gam_frmt = 'milc', 'gwc ', 
addpath /Users/luigi/Work/lib/Matlab

D=length(L);
V=prod(L);
sV=prod(L(1:D-1));
NS=2^(D/2);
T=L(D);
NC=3;
CS=NC*NS;
BB=5;
WITH_TENSOR=0;

%%%%%%%%%%%%%% gamma matrices
[ga gi gv] = gammat(gam_frmt);

%%%%%%%%%%%%%% compute Momenta
[P,Pb,Ph,P2,Pb2,Ph2]= computeMom(L);

%%%%%%%%%%%%%%% initialize 

nconf=size(prop,6);
maxmom=size(prop,5);

prop=permute(prop,[1 3 2 4 5 6]); % c_sink s_sink c_src s_src
prop=reshape(prop,CS,CS,maxmom,nconf);
G=zeros(CS,CS,CS,CS,2*BB);
GF=zeros(CS,CS,CS,CS,2*BB);
avProp=zeros(CS,CS,maxmom);
avProp=mean(prop,4);

V=zeros(CS,CS,4);
A=zeros(CS,CS,4);
GG=zeros(CS,CS,CS,CS,4*BB);
cprop_=zeros(CS,CS);
prop_=zeros(CS,CS);

Gamma=zeros(BB,BB,maxmom,4,nconf);

%%%%%%%%%% Projectors
for mu=1:D
  for  nu=1:D
    sigma(:,:,mu,nu)=0.5*( ga(:,:,mu)*ga(:,:,nu) - ga(:,:,nu)*ga(:,:,mu) );
    sigmat(:,:,mu,nu)=ga(:,:,5)*sigma(:,:,mu,nu);
  end
end
P_AA=zeros(CS,CS,CS,CS);P_VV=zeros(CS,CS,CS,CS);
P_AV=zeros(CS,CS,CS,CS);P_VA=zeros(CS,CS,CS,CS);
P_TT=zeros(CS,CS,CS,CS);P_TTt=zeros(CS,CS,CS,CS);
for a=1:CS
  for b=1:CS
    for c=1:CS
      for d=1:CS
	P_SS(a,b,c,d)=ga(a,b,10)*ga(c,d,10);
	P_PP(a,b,c,d)=ga(a,b, 5)*ga(c,d, 5);
	P_SP(a,b,c,d)=ga(a,b,10)*ga(c,d, 5);
	P_PS(a,b,c,d)=ga(a,b, 5)*ga(c,d,10);
	for mu=1:D
	  P_VV(a,b,c,d)=P_VV(a,b,c,d) + ga(a,b,mu)*ga(c,d,mu);
	  P_AA(a,b,c,d)=P_AA(a,b,c,d) + ga(a,b,mu+5)*ga(c,d,mu+5);
	  P_VA(a,b,c,d)=P_VA(a,b,c,d) + ga(a,b,mu)*ga(c,d,mu+5);
	  P_AV(a,b,c,d)=P_AV(a,b,c,d) + ga(a,b,mu+5)*ga(c,d,mu);
	  for nu=mu+1:D
	    P_TT(a,b,c,d)=P_TT(a,b,c,d) + sigma(a,b,mu,nu)*sigma(c,d,mu,nu);
	    P_TTt(a,b,c,d)=P_TTt(a,b,c,d) + sigma(a,b,mu,nu)*sigmat(c,d,mu,nu);
	  end
	end
      end
    end
  end
end
P_pc_p(:,:,:,:,1)= 1/(64*NC*(NC+1)) * (P_VV + P_AA);
P_pc_m(:,:,:,:,1)= 1/(64*NC*(NC-1)) * (P_VV + P_AA);
P_pc_p(:,:,:,:,2)= 1/(64*(NC^2 -1)) * (P_VV - P_AA) + 1/(32*NC*(NC^2 -1)) * (P_SS - P_PP);
P_pc_m(:,:,:,:,2)= 1/(64*(NC^2 -1)) * (P_VV - P_AA) - 1/(32*NC*(NC^2 -1)) * (P_SS - P_PP);
P_pc_p(:,:,:,:,3)= 1/(32*NC*(NC^2-1)) * (P_VV - P_AA) + 1/(16*(NC^2-1)) * (P_SS - P_PP);
P_pc_m(:,:,:,:,3)=-1/(32*NC*(NC^2-1)) * (P_VV - P_AA) + 1/(16*(NC^2-1)) * (P_SS - P_PP);
P_pc_p(:,:,:,:,4)= (2*NC+1)/(32*NC*(NC^2-1)) * (P_SS + P_PP) - 1/(32*NC*(NC^2-1)) * P_TT;
P_pc_m(:,:,:,:,4)= (2*NC-1)/(32*NC*(NC^2-1)) * (P_SS + P_PP) + 1/(32*NC*(NC^2-1)) * P_TT;
P_pc_p(:,:,:,:,5)=-1/(32*NC*(NC^2-1)) * (P_SS + P_PP) + (2*NC-1)/(96*NC*(NC^2-1)) * P_TT;
P_pc_m(:,:,:,:,5)= 1/(32*NC*(NC^2-1)) * (P_SS + P_PP) + (2*NC+1)/(96*NC*(NC^2-1)) * P_TT;

P_pv_p(:,:,:,:,1)=-1/(64*NC*(NC+1)) * (P_VA + P_AV);
P_pv_m(:,:,:,:,1)=-1/(64*NC*(NC-1)) * (P_VA + P_AV);
P_pv_p(:,:,:,:,2)=-1/(64*(NC^2 -1)) * (P_VA - P_AV) - 1/(32*NC*(NC^2 -1)) * (P_SP - P_PS);
P_pv_m(:,:,:,:,2)=-1/(64*(NC^2 -1)) * (P_VA - P_AV) + 1/(32*NC*(NC^2 -1)) * (P_SP - P_PS);
P_pv_p(:,:,:,:,3)=-1/(32*NC*(NC^2-1)) * (P_VA - P_AV) - 1/(16*(NC^2-1)) * (P_SP - P_PS);
P_pv_m(:,:,:,:,3)= 1/(32*NC*(NC^2-1)) * (P_VA - P_AV) - 1/(16*(NC^2-1)) * (P_SP - P_PS);
P_pv_p(:,:,:,:,4)= (2*NC+1)/(32*NC*(NC^2-1)) * (P_SP + P_PS) - 1/(32*NC*(NC^2-1)) * P_TTt;
P_pv_m(:,:,:,:,4)= (2*NC-1)/(32*NC*(NC^2-1)) * (P_SP + P_PS) + 1/(32*NC*(NC^2-1)) * P_TTt;
P_pv_p(:,:,:,:,5)=-1/(32*NC*(NC^2-1)) * (P_SP + P_PS) + (2*NC-1)/(96*NC*(NC^2-1)) * P_TTt;
P_pv_m(:,:,:,:,5)= 1/(32*NC*(NC^2-1)) * (P_SP + P_PS) + (2*NC+1)/(96*NC*(NC^2-1)) * P_TTt;


%%%%% loop over momenta
for ip=1:maxmom
  ip
  %% S(p) and S^-1(p)
  invProp(:,:) =  inv(avProp(:,:,ip)); %% NSPT inversion

  %%% Zpsi
  uu=zeros(CS,CS);
  for mu=1:4
    uu=uu+(Pb(mu,ip) * (ga(:,:,mu) * invProp(:,:)));
  end
  Zpsi(ip) = -i/(CS) * trace(uu) /Pb2(ip);


  %%%%%%%%%% G_{VV+AA}^{+/-}
  for n=1:nconf
    cprop_(:,:) = ga(:,:,5) * ctranspose(prop(:,:,ip,n)) * ga(:,:,5);
    prop_(:,:)  = prop(:,:,ip,n);
    for mu=1:D
      V_(:,:,mu)=(prop_(:,:) * ga(:,:,mu) * cprop_(:,:)); %% NSPT mult_pert
      A_(:,:,mu)=(prop_(:,:) * ga(:,:,mu) * ga(:,:,5) * cprop_(:,:));
    end
    S_(:,:)=(prop_(:,:) * cprop_(:,:));
    P_(:,:)=(prop_(:,:) * ga(:,:,5) * cprop_(:,:)); 
    if(WITH_TENSOR)
      for mu=1:D
	for nu=1:D
	  T_( :,:,mu,nu)= (prop_(:,:) * sigma( :,:,mu,nu) * cprop_(:,:));
	  Tt_(:,:,mu,nu)= (prop_(:,:) * sigmat(:,:,mu,nu) * cprop_(:,:));
	end
      end
    end

%%%%%%%%%%%%%% G(a,b,c,d,AB) = sum_mu A(a,b,mu) * B(c,d,mu) 
%%%%%%%%%%%%%% GF(a,b,c,d,AB) = sum_mu A(a,c,mu) * B(b,d,mu) 
    ZZ=zeros(CS,CS,CS,CS);
    for mu=1:D
      [XX YY] = ndgrid(V_(:,:,mu),V_(:,:,mu)); %% NSPT mult_pert
      ZZ = ZZ + reshape(XX.*YY,CS,CS,CS,CS);
    end
    G(:,:,:,:,1) = ZZ;
    ZZ=zeros(CS,CS,CS,CS);
    for mu=1:D
      [XX YY] = ndgrid(A_(:,:,mu),A_(:,:,mu));
      ZZ = ZZ + reshape(XX.*YY,CS,CS,CS,CS);
    end
    G(:,:,:,:,2) = ZZ;
    ZZ=zeros(CS,CS,CS,CS);
    for mu=1:D
      [XX YY] = ndgrid(V_(:,:,mu),A_(:,:,mu));
      ZZ = ZZ + reshape(XX.*YY,CS,CS,CS,CS);
    end
    G(:,:,:,:,3) = ZZ;
    ZZ=zeros(CS,CS,CS,CS);
    for mu=1:D
      [XX YY] = ndgrid(A_(:,:,mu),V_(:,:,mu));
      ZZ = ZZ + reshape(XX.*YY,CS,CS,CS,CS);
    end
    G(:,:,:,:,4) = ZZ;

    [XX YY] = ndgrid(S_(:,:),S_(:,:));
    G(:,:,:,:,5) = reshape(XX.*YY,CS,CS,CS,CS);
    [XX YY] = ndgrid(P_(:,:),P_(:,:));
    G(:,:,:,:,6) = reshape(XX.*YY,CS,CS,CS,CS);
    [XX YY] = ndgrid(S_(:,:),P_(:,:));
    G(:,:,:,:,7) = reshape(XX.*YY,CS,CS,CS,CS);
    [XX YY] = ndgrid(P_(:,:),S_(:,:));
    G(:,:,:,:,8) = reshape(XX.*YY,CS,CS,CS,CS);

    if(WITH_TENSOR)
      ZZ=zeros(CS,CS,CS,CS);
      for mu=1:D
	for nu=mu+1:D
	  [XX YY] = ndgrid(T_(:,:,mu,nu),T_(:,:,mu,nu));
	  ZZ = ZZ + reshape(XX.*YY,CS,CS,CS,CS);
	end
      end
      G(:,:,:,:,9) = ZZ;
      ZZ=zeros(CS,CS,CS,CS);
      for mu=1:D
	for nu=mu+1:D
	  [XX YY] = ndgrid(T_(:,:,mu,nu),Tt_(:,:,mu,nu));
	  ZZ = ZZ + reshape(XX.*YY,CS,CS,CS,CS);
	end
      end
      G(:,:,:,:,10) = ZZ;
    end

    GF = permute(G,[1 4 3 2 5]);
    
    GG(:,:,:,:,1) = (G(:,:,:,:,1)-GF(:,:,:,:,1))+(G(:,:,:,:,2)-GF(:,:,:,:,2));%VV+AA|+
    GG(:,:,:,:,2) = (G(:,:,:,:,1)+GF(:,:,:,:,1))+(G(:,:,:,:,2)+GF(:,:,:,:,2));%VV+AA|-
    GG(:,:,:,:,3) = (G(:,:,:,:,1)-GF(:,:,:,:,1))-(G(:,:,:,:,2)-GF(:,:,:,:,2));%VV-AA|+
    GG(:,:,:,:,4) = (G(:,:,:,:,1)+GF(:,:,:,:,1))-(G(:,:,:,:,2)+GF(:,:,:,:,2));%VV-AA|-
    GG(:,:,:,:,5) = (G(:,:,:,:,5)-GF(:,:,:,:,5))-(G(:,:,:,:,6)-GF(:,:,:,:,6));%SS-PP|+
    GG(:,:,:,:,6) = (G(:,:,:,:,5)+GF(:,:,:,:,5))-(G(:,:,:,:,6)+GF(:,:,:,:,6));%SS-PP|-
    GG(:,:,:,:,7) = (G(:,:,:,:,5)-GF(:,:,:,:,5))+(G(:,:,:,:,6)-GF(:,:,:,:,6));%SS+PP|+
    GG(:,:,:,:,8) = (G(:,:,:,:,5)+GF(:,:,:,:,5))+(G(:,:,:,:,6)+GF(:,:,:,:,6));%SS+PP|-
    if(WITH_TENSOR)
      GG(:,:,:,:,9) = (G(:,:,:,:,9)-GF(:,:,:,:,9));%TT|+
      GG(:,:,:,:,10)= (G(:,:,:,:,9)+GF(:,:,:,:,9));%TT|-
    end
   GG(:,:,:,:,11)= (G(:,:,:,:,3)-GF(:,:,:,:,3))+(G(:,:,:,:,4)-GF(:,:,:,:,4));%VA+AV|+
   GG(:,:,:,:,12)= (G(:,:,:,:,3)+GF(:,:,:,:,3))+(G(:,:,:,:,4)+GF(:,:,:,:,4));%VA+AV|-
   GG(:,:,:,:,13)= (G(:,:,:,:,3)-GF(:,:,:,:,3))-(G(:,:,:,:,4)-GF(:,:,:,:,4));%VA-AV|+
   GG(:,:,:,:,14)= (G(:,:,:,:,3)+GF(:,:,:,:,3))-(G(:,:,:,:,4)+GF(:,:,:,:,4));%VA-AV|-
   GG(:,:,:,:,15)= (G(:,:,:,:,8)-GF(:,:,:,:,8))-(G(:,:,:,:,7)-GF(:,:,:,:,7));%PS-SP|+
   GG(:,:,:,:,16)= (G(:,:,:,:,8)+GF(:,:,:,:,8))-(G(:,:,:,:,7)+GF(:,:,:,:,7));%PS-SP|-
   GG(:,:,:,:,17)= (G(:,:,:,:,8)-GF(:,:,:,:,8))+(G(:,:,:,:,7)-GF(:,:,:,:,7));%PS+SP|+
   GG(:,:,:,:,18)= (G(:,:,:,:,8)+GF(:,:,:,:,8))+(G(:,:,:,:,7)+GF(:,:,:,:,7));%PS+SP|-
   if(WITH_TENSOR)
     GG(:,:,:,:,19)= (G(:,:,:,:,10)-GF(:,:,:,:,10));%TTt|+
     GG(:,:,:,:,20)= (G(:,:,:,:,10)+GF(:,:,:,:,10));%TTt|-
   end
   GG=GG./(2);

  %%%%%% Lambda(a',b',c',d') = S(a',a) S(c',c) GG(a,b,c,d) S(b,b') S(d,d')
  Lambda = reshape( invProp * reshape(GG,CS,CS*CS*CS*4*BB), CS,CS,CS,CS,4*BB); %% NSPT mult_pert
  GG= permute(Lambda,[3 2 1 4 5]);
  Lambda = reshape( invProp * reshape(GG,CS,CS*CS*CS*4*BB), CS,CS,CS,CS,4*BB);
  GG= permute(Lambda,[2 3 1 4 5]);
  Lambda = reshape( transpose(invProp) * reshape(GG,CS,CS*CS*CS*4*BB), CS,CS,CS,CS,4*BB);
  GG= permute(Lambda,[4 1 3 2 5]);
  Lambda = reshape( transpose(invProp) * reshape(GG,CS,CS*CS*CS*4*BB), CS,CS,CS,CS,4*BB);
  Lambda= permute(Lambda,[4 2 3 1 5]);

  %%%%%% Gamma(p,q,ip,1,n) =  \sum_a,b,c,d  P_pc_p(a,b,c,d,p) * La_pc_p(b,a,d,c,q);
  La_pc_p = reshape(permute(Lambda(:,:,:,:,[1:2:10]),[2 1 4 3 5]),CS^4,5);
  La_pc_m = reshape(permute(Lambda(:,:,:,:,[2:2:10]),[2 1 4 3 5]),CS^4,5);
  La_pv_p = reshape(permute(Lambda(:,:,:,:,[11:2:20]),[2 1 4 3 5]),CS^4,5);
  La_pv_m = reshape(permute(Lambda(:,:,:,:,[12:2:20]),[2 1 4 3 5]),CS^4,5);

  P_pc_p = reshape(P_pc_p,CS^4,BB);
  P_pc_m = reshape(P_pc_m,CS^4,BB);
  P_pv_p = reshape(P_pv_p,CS^4,BB);
  P_pv_m = reshape(P_pv_m,CS^4,BB);

%  for p=1:BB
%    for q=1:BB
      p=1;
      q=1;
      Gamma(p,q,ip,1,n) =  P_pc_p(:,p)' * La_pc_p(:,q);
      Gamma(p,q,ip,2,n) =  P_pc_m(:,p)' * La_pc_m(:,q);
      Gamma(p,q,ip,3,n) =  P_pv_p(:,p)' * La_pv_p(:,q);
      Gamma(p,q,ip,4,n) =  P_pv_m(:,p)' * La_pv_m(:,q);
%    end
%  end
  
  end % nconf


%  for p=1:BB
%    for q=1:BB
      p=1;
      q=1;
      for o=1:4
      GAM_a(p,q,ip,o)= mean(squeeze(Gamma(p,q,ip,o,:)));
      GAM_d(p,q,ip,o)= std(squeeze(Gamma(p,q,ip,o,:)));

      end
%    end
%  end

  %%% Zo
  nef=4;
  for o=1:4
    DD(:,:,ip,o)=inv(GAM_a(:,:,ip,o));
    Z4F(:,:,ip,o)= Zpsi(ip)^(nef/2) * DD(:,:,ip,o);
  end

  
  
end % ip
