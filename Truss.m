clc
clearvars
close all
%% Inputs
P_1= 1e3; P_2= 2e3;
E_new=[7.3084e10];
Area=[2.38761e-4];
I=[4.32157e-8];
M_Info=[1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;];
Nodes_Coordinates=[0,0;3,0;6,0;9,0;0,4;3,4;6,4;9,4];
Connection=[1,2;2,3;3,4;5,6;6,7;7,8;1,5;2,6;3,7;4,8;1,6;2,7;4,7];
Constr=[0,0;NaN,NaN;NaN,NaN;NaN,0;NaN,NaN;NaN,NaN;NaN,NaN;NaN,NaN];
N_Load=[5,P_1,-P_1;6,0,-P_1;7,0,-P_2;8,0,-P_1];
M_Load=[1,0,0,0,0,0];
%% Get Dj and RL
Nm=size(M_Info,1);
Nn=size(Nodes_Coordinates,1);
Dj(:,1)=reshape(Constr',[2*Nn,1]);
RL=1-isnan(Dj);
%% Get Local and Global Stiffness Matrices
Sj=sparse(2*Nn,2*Nn);
for i=1:Nm
    N1=Connection(i,1);
    N2=Connection(i,2);
    dx=Nodes_Coordinates(N2,1)-Nodes_Coordinates(N1,1);
    dy=Nodes_Coordinates(N2,2)-Nodes_Coordinates(N1,2);
    Lm=sqrt(dx^2+dy^2);
    Cx=dx/Lm;
    Cy=dy/Lm;
    Em=E_new(M_Info(i,1));
    Ame=Area(M_Info(i,2));
    Dof=[2*N1-1,2*N1,2*N2-1,2*N2];    
    Smd(:,:,i)=Em*Ame/Lm*[Cx^2 Cx*Cy -Cx^2 -Cx*Cy;
                          Cx*Cy Cy^2 -Cx*Cy -Cy^2;
                          -Cx^2 -Cx*Cy Cx^2 Cx*Cy;
                          -Cx*Cy -Cy^2 Cx*Cy Cy^2];
    Sj(Dof,Dof)=Sj(Dof,Dof)+Smd(:,:,i);
    RT(:,:,i)= [Cx Cy 0 0;
               -Cy Cx 0 0;
                0 0 Cx Cy;
                0 0 -Cy Cx];
    L(i)=Lm;
end
%% Get Fixed End Actions and Equivalent loads
AE=zeros(2*Nn,1);
Aml=zeros(Nm,4);
Amlb=Aml;
for i=1:size(M_Load,1)
   id=M_Load(i,1);
   Mm=M_Load(i,2);
   Pxm=M_Load(i,3);
   Pym=M_Load(i,4);
   wxm=M_Load(i,5);
   wym=M_Load(i,6);
   N1=Connection(id,1);
   N2=Connection(id,2);
   Lm=L(id);
   Dof=[2*N1-1,2*N1,2*N2-1,2*N2];
   
   Aml(id,1)=-Pxm/2-wxm*Lm/2;
   Aml(id,2)=Mm/Lm-Pym/2-wym*Lm/2;
   Aml(id,3)=Aml(id,1);
   Aml(id,4)=-Mm/Lm-Pym/2-wym*Lm/2;
   Amlb(id,:)=transpose(RT(:,:,id))*transpose(Aml(id,:));
   AE(Dof,1)=AE(Dof,1)-transpose(Amlb(id,:));    
end
%% Get Combined loads
A=zeros(2*Nn,1);
for i=1:size(N_Load,1)
    id=N_Load(i,1);
    Dof=(2*id-1):2*id;
    A(Dof,1)=N_Load(i,2:3);
    Ac=A+AE;
end
%% Renumbering DOF and Get number of possible DOF
ind=[find(RL==0);find(RL==1)];
Re(ind,1)=(1:2*Nn)';
ndpos=length(find(RL==0));
%% Renumbering and Partitioning
Sjre(Re,Re)=Sj;
S=Sjre(1:ndpos,1:ndpos);
Sdr=Sjre(1:ndpos,ndpos+1:2*Nn);
Srd=transpose(Sdr);
Srr=Sjre(ndpos+1:2*Nn,ndpos+1:2*Nn);
Acre(Re,1)=Ac;
Ad=Acre(1:ndpos);
Arl=-Acre(ndpos+1:2*Nn);
Djre(Re,1)=Dj;
Dr=Djre(ndpos+1:2*Nn);
%% Solution
D=S\(Ad-Sdr*Dr);
Djre=[D;Dr];
Dj=Djre(Re)
Ard=Srd*D+Srr*Dr;
AR=Arl+Ard
%% Plotting Deformation
DSF=100;
figure
hold all
for i=1:Nm
    N1=Connection(i,1);
    N2=Connection(i,2);
    Defx=Dj([2*N1-1,2*N2-1],1);
    Defy=Dj([2*N1,2*N2],1);
    Xvec=Nodes_Coordinates([N1 N2],1);
    Yvec=Nodes_Coordinates([N1 N2],2);
    plot(Xvec,Yvec,'--')
    plot(Xvec+DSF*Defx,Yvec+DSF*Defy,'-')
    title('Deformed Shape')
end
%% Plotting Diagrams
xbar1=[0:0.1:0.5]';
xbar2=xbar1+0.5;
Load=zeros(Nm,5);
id=M_Load(:,1);
Load(id,:)=M_Load(:,2:6);
for i=1:Nm
    N1=Connection(i,1);
    N2=Connection(i,2);
    Dof=[2*N1-1,2*N1,2*N2-1,2*N2];
    Dno=Dj(Dof,1);
    Lm=L(i);
    Am(i,:)=Aml(i,:)+transpose(RT(:,:,i)*Smd(:,:,i)*Dno);
    BM1=-Am(i,2)*xbar1*Lm-0.5*Load(i,5)*xbar1.^2*Lm^2;
    BM2=-Am(i,2)*xbar2*Lm-0.5*Load(i,5)*xbar2.^2*Lm^2-Load(i,3)*(xbar2-0.5)*Lm+Load(i,1);
    SF1=-Am(i,2)-Load(i,5)*xbar1*Lm;
    SF2=-Am(i,2)-Load(i,5)*xbar2*Lm-Load(i,3);
    TF1=-Am(i,1)-Load(i,4)*xbar1*Lm;
    TF2=-Am(i,1)-Load(i,4)*xbar2*Lm-Load(i,2);
    figure
    subplot(3,1,1)
    plot(xbar1,TF1,xbar2,TF2)
    xlabel('x_m')
    title(['Thrust Diagram of member ',num2str(i)])
    subplot(3,1,2)
    plot(xbar1,SF1,xbar2,SF2)
    xlabel('x_m')
    title(['Shear Diagram of member ',num2str(i)])
    subplot(3,1,3)
    plot(xbar1,BM1,xbar2,BM2)
    xlabel('x_m')
    title(['Bending Diagram of member ',num2str(i)])
end
Am

