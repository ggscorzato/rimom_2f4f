function [f df] = mean_over_equalP2(g,dg,p2)
% should not be done twice!!

ep=1e-6;
for k=1:length(p2)
  K=find((p2>p2(k)-ep)&(p2<p2(k)+ep));
  f(k)=sum(g(K)./(dg(K)).^2)/sum(1./dg(K).^2);
  df(k) = 1/sqrt(sum(1./dg(K).^2));
end

f=f';
df=df';