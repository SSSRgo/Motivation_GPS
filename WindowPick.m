function [SpkWindow_Nosepoke] = WindowPick(nosepokeon,cellTS,binwidth,upbound,lowbound)
%WINDOWPICK 此处显示有关此函数的摘要
%   此处显示详细说明

% SpkWindow=linspace(0,upbound-lowbound,binwidth).*0;

ii=1;
for i=1:length(nosepokeon)
    if (nosepokeon(i)+lowbound)<0 || (nosepokeon(i)+upbound)>max(cellTS) || isnan(nosepokeon(i))
continue
    end
edges=linspace(nosepokeon(i)+lowbound,nosepokeon(i)+upbound,binwidth);
[SpkWindow_Nosepoke(ii,:),b]=hist(cellTS(find(cellTS<nosepokeon(i)+upbound&cellTS>nosepokeon(i)+lowbound)),edges);
ii=ii+1;
end





% 
% SpkWindow_Reward
% 
% SpkWindow_Nosepoke = nosepokeon;
% outputArg2 = rewardon;
% 
% 


end
% 
% function [SpkWindow_Nosepoke,SpkWindow_Reward] = WindowPick(nosepokeon,rewardon,cellTS,binwidth,upbound,lowbound)
% %WINDOWPICK 此处显示有关此函数的摘要
% %   此处显示详细说明
% 
% % SpkWindow=linspace(0,upbound-lowbound,binwidth).*0;
% 
% ii=1;
% for i=1:length(nosepokeon)
%     if (nosepokeon(i)+lowbound)<0 || (nosepokeon(i)+upbound)>max(cellTS) || isnan(nosepokeon(i))
% continue
%     end
% edges=linspace(nosepokeon(i)+lowbound,nosepokeon(i)+upbound,binwidth);
% [SpkWindow_Nosepoke(ii,:),b]=hist(cellTS(find(cellTS<nosepokeon(i)+upbound&cellTS>nosepokeon(i)+lowbound)),edges);
% ii=ii+1;
% end
% 
% kk=1;
% for k=1:length(rewardon)
%     if (rewardon(k)+lowbound)<0 || (rewardon(k)+upbound)>max(cellTS) || isnan(rewardon(k))
% continue
%     end
% edges=linspace(rewardon(k)+lowbound,rewardon(k)+upbound,binwidth);
% [SpkWindow_Reward(kk,:),b]=hist(cellTS(find(cellTS<rewardon(k)+upbound&cellTS>rewardon(k)+lowbound)),edges);
% kk=kk+1;
% % SpkWindow = discretize(cellTS,edges);
% end
% 
% 
% 
% % 
% % SpkWindow_Reward
% % 
% % SpkWindow_Nosepoke = nosepokeon;
% % outputArg2 = rewardon;
% % 
% % 
% 
% 
% end

