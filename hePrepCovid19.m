function [pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(data,numInt,R0,del)

%data.NNs - column vector of population
%Possible generalsiation to within-sector heterogeneity - one column per subsector

%% POPULATION PARAMETERS:

%Population Density
[n,na]=size(data.NNs);
ntot=n*na;
NN=sum(data.NNs,2);
NNbar=reshape(data.NNs,ntot,1);
NNrep=repmat(NN,na,1);

%Urban and Rural?
%urbrur=0;%Turn in to @home vs @work7/291
%{
if urbrur==1
    hvec=kron(hvec,ones(n,1));
    muvec=kron(muvec,ones(n,1));
    %Kkron=[1,.05;.05,.75];%Sonoma
    Kkron=[1,.05;.05,.85];%Mendocino
else
    Kkron=1;
end
%}
%{
C=[.3827    0.9115    0.0419;
    1.2062    7.3538    0.1235;
    0.1459    0.4810    0.1860];
%}
%C=eye(na);

%Age-Sector Breakdown
lc=4;
adInd=3;
lx=length(data.NNs)-lc;
%xin=ones(lx,1);%DH
xin=data.Xdata(1:63).^(1./data.alpha);
xin=xin';
D=heMakeDs(NN,xin,data,zeros(1,lx));
Dout=D;

%K=heMakeDs(NN,eye(10));
%K=rand(n);
%K=normr(K);
%D=kron(C,K);

pr=struct;
pr.a=data.alpha;%0.411475;

%% DISEASE PARAMETERS:

%from Knock et al. (2020) unless indicated otherwise

%Transmission
pr.R0=R0;%2.6882;%from fitting
pr.red=0.58;%from Byambasuren et al. (2020)

%Latency and Onset
Text=4.6;
pr.sigma=1/Text;
%Tonset=1;
%pr.omega=1/Tonset;

%Case Pathways
pr.p1=0.6;
[ph,pd,~,~]=heParamsAge(data,pr.p1);
% pili=[0.4670,0.4786,0.6590,0.7252]';
% pili=0.505*pili;
% pili=[repmat(pili(adInd),lx,1);pili];
% pr.p2=pili;
% ph=[0.0500    0.0022    0.0350    0.5235]';
ph=[repmat(ph(adInd),lx,1);ph];
% pd=[0.0103    0.0078    0.0361    0.1555]';
pd=[repmat(pd(adInd),lx,1);pd];
% pdeath=0.39;
% pd=pd/sum(pd);
% pd=pdeath*pd;%new data: 39% of hospital cases result in death

%Hospitalisation and Death Rates (proportion and rate combined)
Tsh=4.5;%from Global
%pr.h=ph/Tsh;
Thd=data.Thd;%14; %DH - as Sri Lanka
%pr.mu=pd/Thd;

%Recovery Rates (proportion and rate combined)
Ta=2.1;%asymptomatic
pr.g1=1/Ta;
% pr.gX=1/2.1;%mild
Ts=4;%symptomatic
%pr.g2=(1-ph)/Ts;
Threc=data.Threc;%hospitalised%14;
%pr.g3=(1-pd)/Threc; %DH - as Sri Lanka

%DH:
time_symp=ph.*Tsh+(1-ph).*Ts;
pr.g2=(1-ph)./time_symp;
pr.h=ph./time_symp;
time_hosp=pd.*Thd+(1-pd).*Threc;
pr.g3=(1-pd)./time_hosp;
pr.mu=pd./time_hosp;

%Immunity Loss
Ti=365;%from Global
pr.nu=1/Ti;

%%

Deff=Dout.*repmat(NNbar,1,n*na)./repmat(NNbar',n*na,1);
onesn=ones(ntot,1);

F=zeros(3*ntot,3*ntot);
%F=zeros(4*ntot,4*ntot);
F(1:ntot,ntot+1:end)=[pr.red*Deff,Deff];
%F(1:ntot,ntot+1:end)=[pr.red*Deff,repmat(Deff,1,2)];

vvec=[pr.sigma.*onesn;      pr.g1.*onesn;       (pr.g2+pr.h).*onesn];%g2 and h are vectors
%vvec=[pr.sigma.*onesn;      pr.g1.*onesn;       pr.g2.*onesn;       (pr.gX+pr.h).*onesn];%gX and h are vectors
V=diag(vvec);
V(ntot+1:2*ntot,1:ntot)=    diag(-(1-pr.p1) .*pr.sigma  .*onesn);
V(2*ntot+1:3*ntot,1:ntot)=  diag(-pr.p1     .*pr.sigma  .*onesn);
%V(2*ntot+1:3*ntot,1:ntot)=  diag(-pr.p1     .*(1-pr.p2) .*pr.sigma  .*onesn);
%V(3*ntot+1:4*ntot,1:ntot)=  diag(-pr.p1     .*pr.p2     .*pr.sigma  .*onesn);

GD=F/V;

%Ceff=kron(C,Ckron);%Urb/rural mixing here
%F(1:ntot,1:ntot)=Deff;
%vvec=kron([pr.sigma;pr.g1;pr.g2;pr.gX],ones(ntot,1));
%vvec(end-ntot+1:end)=vvec(end-ntot+1:end)+pr.h;

%{
%HE:
ntot=n*na;
F=zeros(7*ntot,7*ntot);
F(1:ntot,1:ntot)=Deff;%ntot+1:end)=[repmat(2/3*Deff,1,2),repmat(Deff,1,4)];
onesn=ones(ntot,1);
vvec=[pr.sigma*onesn;pr.g1*onesn;pr.omega*onesn;pr.g2*onesn;pr.g2*onesn;pr.h*onesn;pr.h*onesn];
V=diag(vvec);
V(ntot+1:2*ntot,1:ntot)=diag(-pr.sigma*(1-pr.p1)*onesn);
V(2*ntot+1:3*ntot,1:ntot)=diag(-pr.sigma*pr.p1*onesn);
V(3*ntot+1:4*ntot,2*ntot+1:3*ntot)=diag(-(1-pr.p3).*(1-pr.p2));
V(4*ntot+1:5*ntot,2*ntot+1:3*ntot)=diag(-pr.p3.*(1-pr.p2));
V(5*ntot+1:6*ntot,2*ntot+1:3*ntot)=diag(-(1-pr.p4).*pr.p2);
V(6*ntot+1:7*ntot,2*ntot+1:3*ntot)=diag(-pr.p4.*pr.p2);
GD=F/V;
%}

d=eigs(GD,1);%largest in magnitude (+/-) 
R0a=max(d); 
beta=pr.R0/R0a;%beta scales actual R0 to fitted R0

% fun= @(s) s-exp((beta/pr.g2(1))*D*(s-1));
% sinf=fsolve(fun,0.3*ones(size(NN)));
% s=100000*sum(sinf.*NN)/sum(data.NNs);

%% PREPAREDNESS PARAMETERS:

%see country parameter file

%Mitigation Time
% Hm=2*sum(data.Npop)/(10^5);
% ihr=pr.p1*dot(ph,NN)/sum(NN);
% Tg=Text+(1-pr.p1)*Ta+pr.p1*((1-(ihr/pr.p1))*Ts+(ihr/pr.p1)*Tsh);
% r=log(pr.R0)/Tg;
% I0=10^(-3.35);%
% Tm=(log(Hm/(I0*ihr))/r)-60;
% pr.Tm=Tm-0;
%sector closures (x), working-from-home (wfh), NPIs (delta) and self-isolation (p3/p4) implemented

%Adherence to NPIs
pr.betamod=[1,1,del];
% pmod=1.0000;%lockdown delta%from fitting
% pmodIn=1.30*pmod;%post-lockdown delta
% pr.betamod=[1,pmod*ones(1,numInt)];%may be modified in heSingleSim.m or heSwitchSim.m
% pr.betamod=[1,0.5037*ones(1,8),0.8014,0.4633,0.4094,0.4280,0.3971,0.3432,0.2962,0.8584,0.4508,0.45*ones(1,3)];%0.45,0.50,0.60,0.75

%Self-Isolation and Quarantine
pr.p3=0.00;%proportion of asymptomatic self-isolating
pr.p4=0.00;%proportion of symptomatic self-isolating
%pr.p5=0.00;%proportion of severe self-isolating
%
pr.odds=0;%asymptomatic quarantining rate
pr.q1=0;%mild quarantining rate
pr.q2=0;%severe quarantining rate
% pr.g4=1/(1/pr.g2-1/pr.q1);%asymptomatic/mild quarantining recovery rate
% pr.qnew=0;
% if pr.q2>0
%     pr.g4X=1/(1./pr.gX+1./pr.q2);%Vector
% else
%     pr.g4X=pr.q2*0;
% end

%Hospital Capacity
pr.Hmax=data.Hmax;%2017;
%pr.thu=0.90*pr.Hmax;
%pr.thl=0.30*pr.Hmax;
% pr.mu_oc=1.19*pr.mu;%Wilde et al. (2021)
% pr.g3_oc=(1-1.19*pd)/Threc;
%see over-capacity death and recovery rates for vaccinated individuals below

pdoc=heParamsAge_oc(data,pr.p1);
pdoc=[repmat(pdoc(adInd),lx,1);pdoc];
pr.mu_oc=pdoc/Thd;
pr.g3_oc=(1-pdoc)/Threc;

%% VACCINATION PARAMETERS:

vx=struct;%from RTM unless indicated otherwise

%Vaccine 1 (two doses)
vx.hrv1=    1/28;                       %time to develop v-acquired immunity (AstraZeneca)
vx.scv1=    0.58;                       %infection-blocking efficacy
vx.p1v1=    0;                          %disease-blocking efficacy          
vx.hv1=     1-((1-0.90)/(1-vx.scv1));   %severe-disease-blocking efficacy
vx.dv1=     0;                          %death-blocking efficacy
vx.trv1=    0.40;                       %transmission-blocking efficacy
vx.nuv1=    1/Inf;                      %duration of v-acquired immunity

% %Vaccine 1 (one dose)
% vx.hrv1=    1/21;                       %time to develop v-acquired immunity (AstraZeneca)
% vx.scv1=    0.33;                       %infection-blocking efficacy
% vx.p1v1=    0;                          %disease-blocking efficacy          
% vx.hv1=     1-((1-0.75)/(1-vx.scv1));   %severe-disease-blocking efficacy
% vx.dv1=     0;                          %death-blocking efficacy
% vx.trv1=    0.40;                       %transmission-blocking efficacy
% vx.nuv1=    1/Inf;                      %duration of v-acquired immunity

vx.h_v1=    (1-vx.hv1)*ph/Tsh;
vx.g2_v1=   (1-(1-vx.hv1)*ph)/Ts;
vx.mu_v1=   (1-vx.dv1)*pd/Thd;
vx.g3_v1=   (1-(1-vx.dv1)*pd)/Threc;

% vx.mu_ocv1= 1.19*vx.mu_v1;%Wilde et al. (2021)
% vx.g3_ocv1= (1-1.19*(1-vx.dv1)*pd)/Threc;
vx.mu_ocv1=   (1-vx.dv1)*pdoc/Thd;
vx.g3_ocv1=   (1-(1-vx.dv1)*pdoc)/Threc;

% %Vaccine 2
% vx.hrv2=    1/21;                       %time to develop v-acquired immunity (AstraZeneca)
% vx.scv2=    0.65;                       %infection-blocking efficacy
% vx.p1v2=    0;                          %disease-blocking efficacy          
% vx.hv2=     1-((1-0.80)/(1-vx.scv1));   %severe-disease-blocking efficacy
% vx.dv2=     0;                          %death-blocking efficacy
% vx.trv2=    0;                          %transmission-blocking efficacy
% vx.nuv2=    1/Inf;                      %duration of v-acquired immunity

%Administration Rates
ta=data.atimes;%[391,523,548];
arate=1*data.arates;%[9200,16666,112254];

%Population Uptake
Npop=   data.Npop;
NNage=  Npop';%DH %[Npop(1),sum(Npop(2:4)),sum(Npop(5:13)),sum(Npop(14:end))];
pr.up=  data.uptake;%[0.00,0.52*0.90,0.90,0.90];%12:0.52, 14:0.39, 16:0.26, 18:0.13
NNwv=   NNage.*pr.up;

%Rollout
tg=     [];
tg(1)=  ta(1);
rp=     zeros(5,4);
rp(1,:)=[0,0,0,arate(1)];

for i=1:3
    k=5-i;
    j=length(tg);
        
    if tg(j) < ta(2) && (tg(j)+NNwv(k)/arate(1)) < ta(2) %no change in rate
        tg(j+1)=        tg(j)+NNwv(k)/arate(1);
        rp(j+1,k-1)=    arate(1);

    elseif tg(j) < ta(2) && (ta(2)+(NNwv(k)-arate(1)*(ta(2)-tg(j)))/arate(2)) < ta(3) %first change in rate only
        tg(j+1)=        ta(2);
        tg(j+2)=        ta(2)+(NNwv(k)-arate(1)*(tg(j+1)-tg(j)))/arate(2);
        rp(j+1,k)=      arate(2);
        rp(j+2,k-1)=    arate(2);

    elseif tg(j) < ta(2) && j == i %first and second changes in rate
        tg(j+1)=        ta(2);
        tg(j+2)=        ta(3);
        tg(j+3)=        ta(3)+(NNwv(k)-arate(1)*(tg(j+1)-tg(j))-arate(2)*(tg(j+2)-tg(j+1)))/arate(3);
        rp(j+1,k)=      arate(2);   
        rp(j+2,k)=      arate(3);
        rp(j+3,k-1)=    arate(3);

    elseif tg(j) > ta(2) && tg(j) < ta(3) && (tg(j)+NNwv(k)/arate(2)) < ta(3) %no change in rate
        tg(j+1)=        tg(j)+NNwv(k)/arate(2);
        rp(j+1,k-1)=    arate(2);
        
    elseif tg(j) > ta(2) && tg(j) < ta(3) && j == i+1 %second change in rate only
        tg(j+1)=        ta(3);
        tg(j+2)=        ta(3)+(NNwv(k)-arate(2)*(tg(j+1)-tg(j)))/arate(3);
        rp(j+1,k)=      arate(3);
        rp(j+2,k-1)=    arate(3);

    else %no further changes in rate
        tg(j+1)=        tg(j)+NNwv(k)/arate(3);
        rp(j+1,k-1)=    arate(3);

    end    
end

tg_inf=Inf(6,1);%padding
tg_inf(1:length(tg))=tg;
tg=tg_inf;
rp(:,1)=0;%pre-school age group not vaccinated

vx.startp1= tg(1);
vx.startp2= tg(2);
vx.startp3= tg(3);
vx.startp4= tg(4);
vx.startp5= tg(5);
vx.end=     tg(6);

vx.aratep1= rp(1,:)';
vx.aratep2= rp(2,:)';
vx.aratep3= rp(3,:)';
vx.aratep4= rp(4,:)';
vx.aratep5= rp(5,:)';

end

%%

function [phgs,pdgh,Threc,Thd]=heParamsAge(datax,ps)
%%

%nn=[3463,3726,3538,3260,3693,4022,4011,3926,3586,3919,4129,3890,3308,2982,2960,2069,1531+933+414+117+12];%England and Wales, mid-2019 estimate (ONS)
nn=datax.Npop';
%{
nn=[nn(1:16),sum(nn(17:end))];%last age range in Knock et al. (2020) is 80+
ranges=[1,3,9,4];
nntot=[nn(1),sum(nn(2:4)),sum(nn(5:13)),sum(nn(14:end))];
nntot=repelem(nntot,ranges);
%}
ranges=[1,3,9,4];
nntot=repelem(nn,ranges);%DH
nnprop=repelem(1./ranges,ranges);%DH

subs=1:4;
subs=repelem(subs,ranges);

%%

% ps=     0.6;
ihr=    [0.030000 0.002600 0.000840	0.000420 0.000800 0.002600 0.004000	0.006300 0.012000 0.019000 0.023000	0.040000 0.096000 0.100000 0.240000	0.500000 0.500000];
phgs=   ihr./ps;

% picu=   [0.2430 0.2890 0.3380 0.3890 0.4430 0.5030 0.5700 0.6530 0.7560 0.8660 0.9540 1.0000 0.9720 0.8540 0.6450 0.4020 0.1070];
% pg=     1-picu;
% pdicu=  [0.2820 0.2860 0.2910 0.2990 0.3100 0.3280 0.3530 0.3900 0.4460 0.5200 0.6040 0.7050 0.8060 0.8990 0.9690 1.0000 0.9180];
% psd=    1-pdicu;

ifr=    [0.000310 0.000030 0.000010 0.000000 0.000010 0.000040 0.000060 0.000130 0.000310 0.000700 0.001160 0.002760 0.008670 0.012150 0.035120 0.084300 0.096960];
pdgh=   ifr./ihr;

phgs=   accumarray(subs',phgs.*nnprop);
pdgh=   accumarray(subs',pdgh.*nnprop);

%%

% Tgr=    10.7;
% 
% Ttr=    2.5; 
% Ticur=  15.6; 
% Tsdr=   12.2;
% 
% Tgd=    10.3;
% 
% %Ttr=    2.5; 
% Ticud=  11.8;
% 
% %Ttr=    2.5; 
% Ticusd= 7; 
% Tsdd=   8.1;
% 
% Threc=  Tgr*(pg*nn'/sum(nn))+(Ttr+Ticur+Tsdr)*(picu*nn'/sum(nn));
% Thd=    Tgd*(pg*nn'/sum(nn))+(Ttr+Ticud)*((picu.*pdicu)*nn'/sum(nn))+(Ttr+Ticusd+Tsdd)*((picu.*psd)*nn'/sum(nn));

Threc=  NaN;
Thd=    NaN;

end

function pdgh=heParamsAge_oc(datax,ps)
%%

%nn=[3463,3726,3538,3260,3693,4022,4011,3926,3586,3919,4129,3890,3308,2982,2960,2069,1531+933+414+117+12];%England and Wales, mid-2019 estimate (ONS)
nn=datax.Npop';
nn=datax.Npop';
%{
nn=[nn(1:16),sum(nn(17:end))];%last age range in Knock et al. (2020) is 80+
ranges=[1,3,9,4];
nntot=[nn(1),sum(nn(2:4)),sum(nn(5:13)),sum(nn(14:end))];
nntot=repelem(nntot,ranges);
%}
ranges=[1,3,9,4];
nntot=repelem(nn,ranges);%DH
nnprop=repelem(1./ranges,ranges);%DH

subs=1:4;
subs=repelem(subs,ranges);

%%

% ps=     0.6;
ihr=    [0.030000 0.002600 0.000840	0.000420 0.000800 0.002600 0.004000	0.006300 0.012000 0.019000 0.023000	0.040000 0.096000 0.100000 0.240000	0.500000 0.500000];

%ifr=    4*[0.000310 0.000030 0.000010 0.000000 0.000010 0.000040 0.000060 0.000130 0.000310 0.000700 0.001160 0.002760 0.008670 0.012150 0.035120 0.084300 0.096960];
ifr=    [0.0003 0.0003 0.0002 0.0002 0.0003 0.0003 0.0009 0.0009 0.0034 0.0034 0.0073 0.0073 0.0253 0.0253 0.064 0.064 0.1325];
pdgh=   ifr./ihr;

pdgh=   accumarray(subs',pdgh.*nnprop);

end