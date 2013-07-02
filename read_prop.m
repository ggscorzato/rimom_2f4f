function prop = read_prop(confs,mass,L,NC,name_frmt,head_frmt,XP)
% name_frmt= 1: milc, 2: gwc 0:unsplitted_gwc
% head_frmt= 'ildg' or 'gwc '
% XP= 0: config-space, 1: mom-space
% to get all configurations use 
% >./get_conf_num mass 
% and then:
% load confnum; confnum=sort(confnum); confs=confnum([find(diff(confnum));end])';
addpath /Users/luigi/Work/lib/Matlab
D=length(L);
V=prod(L);
sV=prod(L(1:D-1));
NS=2^(D/2);
T=L(D);
RI=2;
B32=4;   % 1 Float = 4 Bytes.
B64=8;
shouldbe=prod(L)*NS*NC*RI; % Expressed in Float32 (=32 bit) or Double.
%skip = 1145 + 127;
skip=33;
checksize=1;
endian='ieee-be';
arch='/Users/luigi/DATA/';

if(XP==0)
%  prop=zeros(NC,NC,NS,NS,sV,T,length(confs));
elseif(XP==1)
  pif=floor([0, T/16,0, L(1)/12, T/12,T/8, L(1)/12,L(1)/8])+1;
  numip=length([pif(1):pif(2)]) * length([pif(3):pif(4)])^3 +...
	length([pif(5):pif(6)]) * length([pif(7):pif(8)])^3;
%  prop=zeros(NC,NC,NS,NS,numip,length(confs));
end


inum=1;
for numc= confs
  if (name_frmt==0)
	filename=strcat(arch,'prop.mass00.',num2str(numc,'%.4d'))
	fid=fopen(filename,'r');
	intro=fread(fid,skip,'uint8=>char');
	ppp=fread(fid,shouldbe*NC*NS,'float64=>float64');
	fclose(fid);
	clear intro
	ppp=reshape(ppp,RI*NC*NS*V,NC,NS);
  end
  for dd=1:NS
    for cc=1:NC
      if (name_frmt==1)
	filename=strcat(arch,'m_',num2str(mass,'%0.3f'),'/config',num2str(numc,'%d'),'_milc_propx.n',...
			num2str(mass,'%0.6f'),'_d',num2str(dd-1,'%d'),'_c',num2str(cc-1,'%d'))
      elseif(name_frmt==2)
	filename=strcat(arch,'prop.mass',num2str(mass,'%2.2d'),'.is',...
			num2str(dd-1,'%d'),'ic',num2str(cc-1,'%d'),'.',num2str(numc,'%.4d'))
%        filename=strcat(arch,num2str(numc,'%.4d'),'/prop.mass',num2str(mass,'%2.2d'),'.is',...
%                        num2str(dd-1,'%d'),'ic',num2str(cc-1,'%d'),'.',num2str(numc,'%.4d'))

      end  
      
      if (checksize==1)
	fid=fopen(filename,'r');
	if(head_frmt=='ildg')
	  intro=fread(fid,3000,'uint8=>char');
	  scidac_point=findstr(transpose(intro),'scidac-binary-data');
	  skip=scidac_point+127; % 127 Byte after scidac-binary-data (found by hand!)
	elseif(head_frmt=='gwc ')
	  [intro countB]=fread(fid,'uint8=>char');
	  skip= countB-shouldbe*B64;
	end
	fclose(fid);
	clear intro
	checksize=0;
	skip
      end
      
      fid=fopen(filename,'r');
      intro=fread(fid,skip,'uint8=>char');
      if(head_frmt=='ildg')
	pp=fread(fid,shouldbe,'float32=>float32');
      elseif(head_frmt=='gwc ')
	pp=fread(fid,shouldbe,'float64',0,endian);
      end
      fclose(fid);
      clear intro

      if (name_frmt==0)
	pp=ppp(:,cc,dd);
      end

      
      if(head_frmt=='ildg')
	pp=reshape(pp,RI,NC,NS,L(1),L(2),L(3),T);
      elseif(head_frmt=='gwc ')
	pp=reshape(pp,RI,NC,NS,T,L(3),L(2),L(1));
	pp=permute(pp,[1 2  3  7 6    5    4]);
      end
      cpp(:,:,:,:,:,:)=squeeze(pp(1,:,:,:,:,:,:) + i * pp(2,:,:,:,:,:,:));
      clear pp

      if(XP==0) 
	prop(:,:,:,:,cc,dd,inum)=reshape(cpp,NC,NS,sV,T);
      elseif(XP==1) % FFT if wanted
	fpp=fft(fft(fft(fft(cpp(:,:,:,:,:,:),[],3),[],4),[],5),[],6);

	ip=1;
	iPTi=pif(1);iPTf=pif(2);iPSi=pif(3);iPSf=pif(4);
	for iPT=iPTi:iPTf
	  for iPY=iPSi:iPSf
	    for iPZ=iPSi:iPSf
	      for iPX=iPSi:iPSf
		prop(:,:,cc,dd,ip,inum) = fpp(:,:,iPX,iPY,iPZ,iPT);
		ip = ip + 1;
	      end
	    end
	  end
	end
	
	iPTi=pif(5);iPTf=pif(6);iPSi=pif(7);iPSf=pif(8);
	for iPT=iPTi:iPTf
	  for iPY=iPSi:iPSf
	    for iPZ=iPSi:iPSf
	      for iPX=iPSi:iPSf
		prop(:,:,cc,dd,ip,inum) = fpp(:,:,iPX,iPY,iPZ,iPT);
		ip = ip + 1;
	      end
	    end
	  end
	end
	clear fpp
      end
    end
  end
  inum = inum + 1;
end

if (XP ==0)
  % convenient rearrangement ? c1 c2 s1 s2 sVol time confnum (c1s1: sink, c2s2: source)
  prop = permute(prop,[1 5 2 6 3 4 7]);
elseif(XP==1)
  % convenient rearrangement ? c1 c2 s1 s2 mom confnum (c1s1: sink, c2s2: source)
  prop = permute(prop,[1 3 2 4 5 6]);
end


% conversion of gamma matrices if needed
if(head_frmt=='ildg')
  [ga gi gv tgi tgv] = gammat('gwc ');
  if (XP==0)

    for inum=1:length(confs)
      for t=[1:T]  %% SEGNO
	for s1=1:NS
	  for s2=1:NS
	    temp(:,:,s1,s2,:)=gv(2,s1)*prop(:,:,gi(2,s1),tgi(2,s2),:,t,inum)*tgv(2,s2);
	  end
	end
	prop(:,:,:,:,:,t,inum)=temp(:,:,:,:,:);
      end
    end
  elseif(XP==1)
    for inum=1:length(confs)
      for s1=1:NS
	for s2=1:NS
	  temp(:,:,s1,s2,:)=gv(2,s1)*prop(:,:,gi(2,s1),tgi(2,s2),:,inum)*tgv(2,s2);
	end
      end
      prop(:,:,:,:,:,inum)=temp(:,:,:,:,:);
    end
  end
  
end
