masses=[4];
allconfs=[1400 1450 1500 1550 1600 1650 1700 1750 1800 1850 1900 1950 2000 2050 2100 2150 2200 2250 2300 2350 ...
2400 2450 2500 2550 2600 2650 2700 2750 2800 2850 2900 2950 4150 4200 4250 4300 4350 4400 4450 4500 4550 4600 ...
4650 4700 4750 4800 4850 4900 4950 5000 5050 5100 5150 5200];

%eval(['load results.mat']);

NC=3;
L=[24 24 24 48];
name_frmt=2;
head_frmt='gwc ';
gam_frmt='gwc ';
arch='/Users/luigi/DATA/momprops'
%arch='/afs/ifh.de/group/nic/scratch/pool1/scorzato/momprops';
XP=1;
[P,Pb,Ph,P2,Pb2,Ph2]= computeMom(L);

for mm=1:length(masses)
  mm
  if (masses(mm)==9)
    confs=allconfs([1:end-2, end]);
  else
    confs=allconfs;
  end
  for cc=1:length(confs)
    file=strcat(arch,'/momprop_',num2str(confs(cc),'%.4d'),'.',num2str(masses(mm),'%.2d'),'.mat');
    eval(['load ',file]);
    fprop(:,:,:,:,:,cc)=prop;
  end
     
  clear prop
%  prop=read_prop(confs,mm,L,NC,name_frmt,head_frmt,XP);
%  [avProp(:,:,:,mm),Gamma2(:,:,mm),Zpsi(:,mm),Zo(:,:,mm)] = computeZ(fprop,L,gam_frmt);
  [Zpsi1_b,Z4F_b(:,:,:,:,mm),Gamma4_b(:,:,:,:,mm),dGamma4_b(:,:,:,:,mm)] = ...
      computeZ4F_onlyBK_we_bak(fprop,L,gam_frmt);
  [Zpsi1,Z4F(:,:,:,:,mm),Gamma4(:,:,:,:,mm),dGamma4(:,:,:,:,mm)] = computeZ4F_onlyBK_we(fprop,L,gam_frmt);
  
  clear fprop
  eval(['save results1.mat']);

end

