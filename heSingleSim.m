function [g,h]=heSingleSim(xoptim)

    load('SR35.mat','data');
    numInt=length(xoptim)/length(data.G);

    [pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(data,numInt,2.7947,[0.5241*ones(1,8),0.9009,0.4402,0.4197,0.4100,0.4084,0.3277,0.2920,0.8650,0.4510,0.4510,...
                                                                                 0.4510*ones(1,3),1.5000*ones(1,3)]);%0.4510,0.6039,0.7842,1.2748,1.5000
    
    %%
    
    % pr.sw=0;%switching off
    
    tvec=[-41.1072,1,32,61,92,122,153,183,214,245,275,306,336,367,398,426,457,487,518,548,579,610,640,671,701,731];           
    data.wfhAv=[zeros(3,length(data.G));repmat(data.wfhAv,numInt-2,1)];
    woptim=xoptim.^(1/pr.a);
    %woptim(35*19+32)=0.1;
    %%
    
    [f,g,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,woptim,tvec,0,data);

    plotMultiOutd_sr(f,woptim,tvec,data,pr);
    
    %%
    ind=find(f(:,1)>tvec(end-6),1)-1;
    h(1)=f(end,5)-f(ind,5);%deaths
    h(1)=round(h(1));
    %h(1)=h(1)/(sum(data.Npop)/(10^5));%deaths (per 100k)
    
%     fullx=[ones(length(data.G),1),reshape(xoptim,length(data.G),numInt)];
%     
%     dgva=repmat((12/365)*data.obj,1,length(tvec)-1);
%     
%     h(2)=(tvec(end)-tvec(1))*sum((12/365)*data.obj)-(tvec(2:end)-tvec(1:end-1))*sum(fullx.*dgva,1)';%GDP loss ($, million)
%     %h(2)=h(2)/1000;%GDP loss ($, billion)%h(2)=round(100*h(2)/2/(12*sum(data.obj)),4);%GDP loss (% per annum)

    fullx=reshape(xoptim(35*18+1:end),length(data.G),6);
    
    mgva=repmat(data.obj,1,6);
    
    h(2)=sum(mgva,'all')-sum(fullx.*mgva,'all')';%GDP loss ($, million)
    h(2)=round(h(2),2);
    
    figHandles=findobj('Type','figure');
    figure(figHandles(2));
    title(['\textbf{Deaths:} ' num2str(h(1)) ' \big/ \textbf{GDP Loss:} \$' num2str(h(2)) ' (million)']);

end