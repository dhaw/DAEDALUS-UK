function f=ukRoadmapMakeXsSimple(Xin)
%Input grid

%X - Daedalus input for roadmap periods
%x - monthly values (from data)
Xout=Xin;
lx=size(Xin,1);%Number of sectors

repCol=[3,1,19,1,2,1,3,5,1,4,3,1,5,4,1,1,2,2,3,1];
repColSum=[1,cumsum(repCol)+1];

A=zeros(6,6);
A(1,1)=28;
A(2,1:3)=[4,21,3];
A(3,3:4)=[11,19];
A(4,4:5)=[16,15];
A(5,5)=30;
A(6,5:6)=[18,13];
a=[28,31,30,31,30,31]';
[a1,a2]=size(A);
A=A./repmat(a(1:a1),1,a2);
undoA=inv(A);

%nullA=null(A);
%undoA=A'\(A'*A);

for i=1:length(repCol)
    %Simpler method:
    %Feb to July or equiv periods
    %
    numSect=repColSum(i+1)-repColSum(i);
    xi=Xin(repColSum(i),13:18)';
    x0=undoA*xi;
    %xi=C*xi;
    Xout(repColSum(i):repColSum(i+1)-1,13:18)=repmat(x0',numSect,1);%B*x0;
    %}
end
f=Xout;
end
