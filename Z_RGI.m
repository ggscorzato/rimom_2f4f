load results2.mat
load results4we.mat

load bare_masses.txt

mu1=2.2; % scale in GeV to compute R^RI(mu)
mui=2.5;  % The scale above which the plateau for Z^RGI and Z^RI-MOM is assumed
mue=3.5; 
nconf=54;
cut=0.4;
obs=[1 3 4]; % S P V A

P2r=P2(find(P2>cut));
I2r=find(P2>cut);

a=(0.086/197.3)*1e3; % inverse GeV
r0_a=5.22;
P2atmu1= (a*mu1)^2;
P2atmui= (a*mui)^2;
P2atmue= (a*mue)^2;
p2_range=find((P2>P2atmu1-0.1) & (P2<P2atmu1+.1));
p2_rgi=find((P2>=P2atmui) & (P2<=P2atmue));
p2_rimom=find((P2>=P2atmui) & (P2<=P2atmue));

bm_done=[3:14]; % index goes from 1 to 21.
bm_chose=[3:17];

for mm = bm_done
fv=real(Gamma2(:,3,mm));
dfv=real(dGamma2(:,3,mm));

go=real(squeeze(Gamma4(1,1,:,1,mm)));
dgo=real(squeeze(dGamma4(1,1,:,1,mm)))./sqrt(nconf); % in future this should be done inside computeZ4, as for Z2
[fo dfo]=mean_over_equalP2(go,dgo,P2); % in future this should be done inside computeZ4, as for Z2


R_RI(:,mm) = fv.^2./fo;
dR_RI(:,mm) = R_RI(:,mm) .* sqrt((2*dfv./fv).^2 + (dfo./fo).^2);

end


%errorbar(P2r,R_RI(I2r,14),dR_RI(I2r,14),'k.')


for ip = 1:length(P2)
  [pol, str]=polyfit(bare_masses(bm_done,2),R_RI(ip,bm_done)',1);
  [R_RI0(ip), err]= polyval(pol,0,str);
  [dR_RI0(ip),which]=max([err,mean(dR_RI(ip,bm_done))]); %  dR_RI0(ip)=err;

  [pol, str]=polyfit(bare_masses(bm_chose,2),Zpsi(ip,bm_chose)',1);
  [Zpsi0(ip), err]= polyval(pol,0,str);
  [dZpsi0(ip),which]=max([err,mean(dZpsi(ip,bm_chose))]);
  for o=obs
    [pol, str]=polyfit(bare_masses(bm_chose,2),real(squeeze(Zo(ip,o,bm_chose))),1);
    [Zo0(ip,o), err]= polyval(pol,0,str);
    [dZo0(ip,o),which]=max([err,mean(dZo(ip,o,bm_chose))]);
  end
  for o=obs
    [pol, str]=polyfit(bare_masses(bm_chose,2),real(squeeze(Gamma2(ip,o,bm_chose))),1);
    [Gamma20(ip,o), err]= polyval(pol,0,str);
    [dGamma20(ip,o),which]=max([err,mean(dGamma2(ip,o,bm_chose))]);
  end
end


mur0= sqrt(P2).*r0_a;
u = URI(mur0,2);

z = R_RI0 .* u;
dz = dR_RI0 .* u;

R_RI_mu=mean(R_RI0(p2_range));
dR_RI_mu=mean(dR_RI0(p2_range));

ZS_RI_mu=mean(Zo0(p2_range));
dZS_RI_mu=mean(dZo0(p2_range));

z_rgi=mean(z(p2_rgi));   % dirty !! weight in the same way different momenta
dz_rgi=mean(dz(p2_rgi));

for ii=1:length(P2)
ZS_rgi(ii,1) = RG_ZPS(sqrt(P2(1,ii)),1,3,0.245,Zo0(ii,1),2);
dZS_rgi=dZo0(:,1);
end

for o=obs
  [polo(:,o), str]=polyfit(P2(p2_rimom),Zo0(p2_rimom,o)',1);
  [z_rimom(o), err]= polyval(polo(:,o),0,str);
  [dz_rimom(o),which]=max([err,mean(dZo0(p2_rimom,o))]);
end
for o=obs
  [polog(:,o), str]=polyfit(P2(p2_rimom),Gamma20(p2_rimom,o)',1);
  [g_rimom(o), err]= polyval(polog(:,o),0,str);
  [dg_rimom(o),which]=max([err,mean(dGamma20(p2_rimom,o))]);
end
[pol_psi, str]=polyfit(P2(p2_rimom),Zpsi0(p2_rimom),1);
[zpsi_rimom, err]= polyval(pol,0,str);
[dzpsi_rimom,which]=max([err,mean(dZpsi0(p2_rimom))]);


figure(1)
hold
errorbar(P2r,R_RI0(I2r),dR_RI0(I2r),'kx')
%for mm= bm_done
%c=[0 1 1/mm];
%errorbar(P2r,R_RI(I2r,mm),dR_RI(I2r,mm),'Color',c,'LineStyle','-')
%end
errorbar(P2r,z(I2r),dz(I2r),'bx')

on1=ones(length(p2_range),1);
on2=ones(length(p2_rgi),1);

%errorbar(P2(p2_range),R_RI_mu*on1,dR_RI_mu*on1,'r')
errorbar(P2(p2_rgi),z_rgi*on2,dz_rgi*on2,'b')

figure(2)
hold
errorbar(P2r,Zpsi0(I2r),dZpsi0(I2r),'.k')
sP2=sort(P2);
plot(sP2,pol_psi(1)*sP2+pol_psi(2),'-k')

figure(3)
hold
errorbar(P2r,Zo0(I2r,o),dZo0(I2r,o),'.k')
plot(sP2,polo(1,o)*sP2+polo(2,o),'-k')
o=4;
errorbar(P2r,Zo0(I2r,o),dZo0(I2r,o),'.r')
plot(sP2,polo(1,o)*sP2+polo(2,o),'-r')

figure(4)
hold
o=1;
errorbar(P2r,Gamma20(I2r,o),dGamma20(I2r,o),'.g')
plot(sP2,polog(1,o)*sP2+polog(2,o),'-g')
o=3;
errorbar(P2r,Gamma20(I2r,o),dGamma20(I2r,o),'.k')
plot(sP2,polog(1,o)*sP2+polog(2,o),'-k')
o=4;
errorbar(P2r,Gamma20(I2r,o),dGamma20(I2r,o),'.r')
plot(sP2,polog(1,o)*sP2+polog(2,o),'-r')

figure(5)
o=1;
hold
errorbar(P2r,Zo0(I2r,o),dZo0(I2r,o),'.b')
%plot(sP2,polo(1,o)*sP2+polo(2,o),'-g')
%errorbar(P2r,ZS_rgi(I2r),dZS_rgi(I2r),'.r')
o=3;
