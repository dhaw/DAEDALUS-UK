function [f,g,tvec,rhohat]=beRunCovid19(pr,be,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,Xit,tvec,plotTau,data)%,pDeff=Dout.*repmat(NNbar,1,n*na)./repmatmodIn)

pr.Xmu=mean(Xit,2);

%% MODEL TYPE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trivialEconomy=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Inputs up to beta are outputs from hePrepCovid19
%Xit - column vector with proportion of each sector open at each
%intervention point. 
%tvec - vector of time points including tvev(1)=t0, tvec(2)=lockdown start,
%tvec(end)=end of simulation. 
%plotTau=1 to plot output. 
%Note: this is set to 0 in any optimisation protocol to avoid a crash due
%to rendering loads of images!

%%

lt=length(tvec);
numseed=7;
phi1=1; phi2=.1;
pr.tlag=32;%1st Feb peak
tauend=1;

%% PANDEMIC PREPAREDNESS:

lc=4;%*age
adInd=3;%Age group of adults *age
lx=size(Xit,1);%length(data.G);%Number of sectors %*1sector

%X's in include pre-lockdown
%BUT making of X's assumes Feb-April economy at lockdown date
%For strict onset of lockdown, rather than averaged-out March
XitMat=reshape(Xit,lx,lt-1);
pr.x=XitMat;
NNvec=repmat(NNbar(1:lx),1,lt-1).*XitMat;%Assumes pre-lockdown=fully open
NNworkSum=sum(NNvec,1);
NNvec(lx+1:lx+lc,1:lt-1)=repmat(NNbar(lx+1:lx+lc),1,lt-1);
NNvec(lx+adInd,:)=sum(NNbar([1:lx,lx+adInd]))-NNworkSum;

%D=Dout;
numComps=length(NN)*be.numGroups;
Dvec=zeros(numComps,numComps,lt-1);
NNvecHold=zeros(size(NNvec));
NNvecHold=repmat(NNvecHold,2,1);
ntot=2*size(be.propsB,1);%xx OVERWRITTEN HERE
pOrder=zeros(ntot,lt-1);
for i=1:lt-1
    Dvec(:,:,i)=pr.betamod(i)*beMakeDs(NNvec(:,i),XitMat(:,i),data,data.wfhAv(i,:),be,be.propsB(:,i));
    pOrderi=reshape([be.propsB(:,i)';1-be.propsB(:,i)'],ntot,1);
    NNvecHold(:,i)=kron(NNvec(:,i),ones(2,1)).*pOrderi;
    pOrder(:,i)=pOrderi;
end
NNvec=NNvecHold;%%xx Not ncessarily integer populations
%%
A1=zeros(n,lt);
A2=A1;
seed=10^-(numseed);%*NNprob;
seedvec=seed*ones(ntot,1); %zeros(ntot,1); seedvec(2*n+1:3*n)=seed*ones(n,1);
%%
%SIMULATE (UP TO ATTACK RATES):
for t=1:tauend%t=tau
[DEout,Rout,tvec,rhohat]=beSimCovid19vax(pr,be,vx,beta,tvec,Dvec,n,ntot,NNvec,phi1,phi2,seedvec,NNvec(:,1),t,plotTau,pOrder);%S0=NNbar (2nd last arg)
%[DEout,Rout]=heSimCovid19vax(pr,beta,tvec,Dvec,n,ntot,NNvec,phi1,phi2,seedvec,NNvec(:,1),t,plotTau);%S0=NNbar (2nd last arg)
%nu=Zsol-Z0;
A1=DEout;%(:,t)=sum(reshape(Zsol./NN0,n,na),2);%Prop immune for spatial cell (before antigenic drift)
A2=Rout;%(:,t)=sum(reshape(nu./NN0,n,na),2);%AR for spatial cell
end
%%
%}
f=A1;
g=A2;
end