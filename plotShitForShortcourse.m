function f=plotShitForShortcourse
x=(-3:1:3);
y1=6*x-2;
y2=1-x/2;
y3=-x-1;
fs=12; lw=2;
cmap=lines(7);
maxy=max([y1,y2,y3]);
miny=min([y1,y2,y3]);
figure
hold on
plot([0,0],[miny,maxy],'k-','linewidth',1)
plot([x(1),x(end)],[0,0],'k-','linewidth',1)
h1=plot(x,y1,'-','linewidth',lw,'color',cmap(1,:));
h2=plot(x,y2,'-','linewidth',lw,'color',cmap(2,:));
h3=plot(x,y3,'-','linewidth',lw,'color',cmap(3,:));
legend([h1,h2,h3],'Q1','Q2','Q3')
set(gca,'fontsize',fs)
xlabel('x')
ylabel('y','rot',0)
axis tight
box on
grid on


x=(-1:.01:1);
y1=.5.^x;
y2=8.^x;
y3=2.^(-3*x)+1;
maxy=max([y1,y2,y3]);
miny=min([y1,y2,y3]);
figure
hold on
plot([0,0],[0,maxy],'k-','linewidth',1)
plot([x(1),x(end)],[0,0],'k-','linewidth',1)
hold on
h1=plot(x,y1,'-','linewidth',lw,'color',cmap(4,:));
h2=plot(x,y2,'-','linewidth',lw,'color',cmap(5,:));
h3=plot(x,y3,'-','linewidth',lw,'color',cmap(6,:));
legend([h1,h2,h3],'Q1','Q2','Q3')
set(gca,'fontsize',fs)
xlabel('x')
ylabel('y','rot',0)
axis tight
box on
grid on