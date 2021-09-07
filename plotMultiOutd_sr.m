function f=plotMultiOutd_sr(f1,x1,tvec,data,pr)

numSectors=length(data.G);
numPeriods=length(tvec)-1;
dodiff=1;

%%

f=figure('Units','centimeters','Position',[0 0 30 10]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',12);fs=12; lw=2;
ax=gca;
ax.Position=[0.05 0.175 0.90 0.75];
cmap=lines(2);
thresh=pr.Hmax;
numThresh=length(thresh);
t1=f1(:,1); 
s1=f1(:,2); 
I1=f1(:,3);
h1=f1(:,4);
d1=f1(:,5);
v1=f1(:,6);%heSimCovid19 output
v2=f1(:,7);%heSimCovid19 output
v3=f1(:,8);%heSimCovid19 output
v4=f1(:,9);%heSimCovid19 output
if dodiff==1
    inc1=-diff(s1,1);%f1=Sout
    tdiff=diff(t1,1);
    inc1=inc1./tdiff;%repmat - older version?
end
hold on

scal1=sum(data.Npop)/(10^5);
scal2=sum(data.Npop)/(10^6);
scal3=sum(data.Npop)/(10^7);
scal4=sum(data.Npop)/(10^8);

%maxY=max([5*max(I1/scal1)/4,5*max(d1)/4,100000,pr.Hmax,5*max(h1)/4]);%inc2;h2
%maxY=max([100000,(5/4)*d1'/scal4,(5/4)*thresh/scal4,(5/4)*h1'/scal4,(5/4)*I1'/scal3]);%inc2;h2
maxY=max([(5/4)*(v1+v2+v3+v4)'/scal1,(5/4)*d1',(5/4)*thresh,(5/4)*h1',(5/4)*I1'/scal3,25000]);

for i=2:length(tvec)-1 
    %for i=2:length(tvec)-1
    plot(tvec(i)*[1,1],[0,maxY],'k-','linewidth',0.01)
end
a=[tvec(end-6) tvec(end-6) tvec(end) tvec(end)];
b=[-100 maxY+100 maxY+100 -100];
fill(a,b,'yellow','linewidth',0.01,'facealpha',.2);
for j=1:numThresh
    plot([0,tvec(end)],[thresh(j),thresh(j)],':','linewidth',lw,'color',.5*[1,1,1])
end

Npop=   data.Npop;
NNage=  [Npop(1),sum(Npop(2:4)),sum(Npop(5:13)),sum(Npop(14:end))];
if ((pr.up(1)-0.02<v1(end)/NNage(1))&&(v1(end)/NNage(1)<pr.up(1)+0.02))&&...
   ((pr.up(2)-0.02<v2(end)/NNage(2))&&(v2(end)/NNage(2)<pr.up(2)+0.02))&&...
   ((pr.up(3)-0.02<v3(end)/NNage(3))&&(v3(end)/NNage(3)<pr.up(3)+0.02))&&...
   ((pr.up(4)-0.02<v4(end)/NNage(4))&&(v4(end)/NNage(4)<pr.up(4)+0.02))
    disp('');
else
    error('Rollout Error!');
end 

hh9=plot(t1,10000*v1/NNage(1),'-','linewidth',lw,'color','green');
hh8=plot(t1,10000*v2/NNage(2),':','linewidth',lw,'color','green');
hh7=plot(t1,10000*v3/NNage(3),'--','linewidth',lw,'color','green');
hh6=plot(t1,10000*v4/NNage(4),'-.','linewidth',lw,'color','green');
hh5=plot(t1,(v1+v2+v3+v4)/scal1,'-','linewidth',lw,'color','green');
hh4=plot(t1,d1,'-','linewidth',lw,'color','black');
hh3=plot(t1,h1,'-','linewidth',lw,'color','magenta');
hh2=plot(t1,I1/scal3,'-','linewidth',lw,'color','red');
%hh1=plot(t1(1:end-1)+0.5*tdiff,inc1,'-','linewidth',lw,'color','yellow');
%hh1=plot(t1,maxY*s1/sum(data.Npop),'-','linewidth',lw,'color','blue');

points=tvec+5;
pointsy=.95*maxY;
txt={'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30',...
     '31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','46','47','48','49','50','51','52','53','54','55','56','57','58','59','60'};
%text(5,pointsy,'PRE','fontsize',fs);
%text(tvec(2)+5,pointsy,'LD','fontsize',15);
for i=10:numPeriods
    text(points(i),pointsy,txt{i-1},'fontsize',fs);
end
axis ([tvec(10),tvec(end),0,maxY])
xlabel('Time','FontSize',fs);
ylabel('Number','FontSize',fs);%yvar
vec_pos=get(get(gca,'ylabel'),'Position');
set(get(gca,'ylabel'),'Position',vec_pos+[-2.5 0 0]);
set(gca,'FontSize',fs);

%xticks([1,32,61,92,122,153,183,214,245,275,306,336,367,398])
%xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb'})
set(gca,'xtick',[[1,32,61,92,122,153,183,214,245,275,306,336],366+[1,32,60,91,121,152,182,213,244,274,305,335],366+365+[1,32,60,91,121,152,182,213,244,274,305,335]]);
set(gca,'xticklabels',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec',...
                       'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec',...
                       'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
xtickangle(45);
ax = gca;
ax.YAxis.Exponent = 3;
%legend([hh1,hh2],'Inc.','Hosp. occ.','location','west')
legend([hh2,hh3,hh4,hh5],'Prevalence (per 10m)','Hospital Occupancy','Deaths','Vaccinated (per 100k)','location','west');%'Position',[-0.29 0.27 1 1]);
%legend([hh3,hh4],'Hospital Occupancy','Deaths','location','west');%'Position',[-0.29 0.27 1 1]);
%legend([hh1,hh2,hh3,hh4],'Inc
%legend([hh1,hh2,hh3,hh4],'Inc. (xmin)','HC (xmin)','Inc. (open)','HC (open)','location','west')
grid on;
%grid minor
box on
hold off

%%

f=figure('Units','centimeters','Position',[0 0 10 10]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'defaultColorbarTickLabelInterpreter','latex');
set(f,'DefaultAxesFontSize',12);
ax=gca;
ax.Position=[0.15 0.10 0.80 0.80];

% if pr.sw==0
    x1=[ones(numSectors,1);x1];
    x=reshape(x1,numSectors,numPeriods);
% else
%     x1=repmat(reshape(x1,numSectors,2),1,length(tvec));
%     x=[ones(numSectors,1),x1(:,1:length(tvec)-2)];
% end

xvec=(1:numPeriods)';xvec=xvec-.5;
xlabels2=({'PRE','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24'});
yvec=(1:5:numSectors)';
ylab=num2str(yvec);
yvec=yvec-.5;
colormap hot;
imagesc(x)
set(gca,'YDir','normal')
xlim([13.5,25.5]);
xlabel('Period')
ylabel('Sector')
vec_pos=get(get(gca,'ylabel'),'Position');
set(get(gca,'ylabel'),'Position',vec_pos+[-0.75 0 0]);
set(gca,'fontsize',fs,'xtick',xvec,'xticklabels',xlabels2,'ytick',yvec,'yticklabels',ylab);%{'PRE','LD',xlab(3:end,:)}
axis square;%([0,numPeriods,.5,63.5])%yvec(1),yvec(end)+1])
xtickangle(45);
caxis([0,1])
colorbar
grid on
box on

end