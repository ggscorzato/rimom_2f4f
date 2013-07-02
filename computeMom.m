function [P,Pb,Ph,P2,Pb2,Ph2]= computeMom(L)
% P(i,n)= 2 pi k_i/L_i; k_i=0:L_i-1
% Pb= sin(P)
% Ph= 2 sin(P/2)
% Pb2 = sum_i Pb(i,:).^2; Ph2 = sum_i Ph(i:).^2; 

PHASES=2.*pi/L(1);
PHASET=2.*pi/L(4);
T=L(4);
%pif=floor([0, T/16,0, L(1)/12, T/12,T/8, L(1)/12,L(1)/8])+1;
pif=[0,3,0,2,4,6,2,3]+1;

ip=0;
%%%%first part
iPTi=pif(1);
iPTf=pif(2);
iPSi=pif(3);
iPSf=pif(4);
for iPT=iPTi:iPTf
  for iPY=iPSi:iPSf
    for iPZ=iPSi:iPSf
      for iPX=iPSi:iPSf
	ip=ip+1;
	P(1,ip)=(PHASES*(iPX-1));
	P(2,ip)=(PHASES*(iPY-1));
	P(3,ip)=(PHASES*(iPZ-1));
	P(4,ip)=(PHASET*(iPT-1));
      end
    end
  end
end
%%%%%%%%%%%%% second part
iPTi=pif(5);
iPTf=pif(6);
iPSi=pif(7);
iPSf=pif(8);
for iPT=iPTi:iPTf
  for iPY=iPSi:iPSf
    for iPZ=iPSi:iPSf
      for iPX=iPSi:iPSf
	ip=ip+1;
	P(1,ip)=(PHASES*(iPX-1));
	P(2,ip)=(PHASES*(iPY-1));
	P(3,ip)=(PHASES*(iPZ-1));
	P(4,ip)=(PHASET*(iPT-1));
      end
    end
  end
end

P2=sum(P.^2,1);
Pb=sin(P);
Ph=2*sin(P/2);
Pb2=sum(Pb.^2,1);
Ph2=sum(Ph.^2,1);
