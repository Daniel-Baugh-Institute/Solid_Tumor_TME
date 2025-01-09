Y_01=lhsdesign(500,10);
Y_02=5000*lhsdesign(500,10);
Y_0=[Y_01;Y_02];
Y_0(:,5)=0;
tspan=[0 70000];
P=ST_parameters_Immune_Desert_Gompertz;
X_ST=[];
K_IL2_TK=[0.25 1.5];
C1=zeros(length(Y_0),3);
for l=1:length(K_IL2_TK)
P(23)=K_IL2_TK(l);
figure
ax1=axes();
for k=1:1:length(Y_0)
[t_pre,x_pre]=ode23s(@(t,y)ST_mod_Immune_Desert_Gompertz(t,y,0,0,0,P),tspan,Y_0(k,:)');
plot(x_pre(:,1)/10000,log2(x_pre(:,4)))
xlabel('PDL1- Tumor cells population','FontSize',20,'FontWeight','normal','FontName','Palatino Linotype')
ylabel('Killer T cells population','FontSize',20,'FontWeight','normal','FontName','Palatino Linotype')
set(ax1,'FontName','Palatino Linotype','FontSize',19,'FontWeight','normal','LineWidth',1);
hold on
X_ST(k,:)=x_pre(length(x_pre(:,1)),:);
scatter(X_ST(k,1)/10000,log2(X_ST(k,4)),800,[0 0 0],"x")
box off
end
end

for l=1:length(K_IL2_TK)
P(23)=K_IL2_TK(l);
figure
ax1=axes();
for k=1:1:length(Y_0)
[t_pre,x_pre]=ode23s(@(t,y)ST_mod_Immune_Desert_Gompertz(t,y,0,0,0,P),tspan,Y_0(k,:)');
plot(log2(x_pre(:,10)),log2(x_pre(:,4)))
xlabel('IL-2 concentration','FontSize',20,'FontWeight','normal','FontName','Palatino Linotype')
ylabel('Killer T cells population','FontSize',20,'FontWeight','normal','FontName','Palatino Linotype')
set(ax1,'FontName','Palatino Linotype','FontSize',19,'FontWeight','normal','LineWidth',1);
hold on
X_ST(k,:)=x_pre(length(x_pre(:,1)),:);
scatter(log2(X_ST(k,10)),log2(X_ST(k,4)),800,[0 0 0],"x")
box off
end
end

X_STT=[];
for l=1:length(K_IL2_TK)
P(23)=K_IL2_TK(l);
ax1=axes();
for k=1:1:length(Y_0)
[t_pre,x_pre]=ode23s(@(t,y)ST_mod_Immune_Desert_Gompertz(t,y,0,0,0,P),tspan,Y_0(k,:)');
X_STT(k,:,l)=x_pre(length(x_pre(:,1)),:);
end
end
L11=find(X_STT(:,4,1)<0.08);
L12=find(X_STT(:,4,1)>0.08);
L21=find(X_STT(:,4,2)<0.2);
L22=find(X_STT(:,4,2)>0.2);
C1=zeros(length(X_STT(:,4,1)),3);
C2=C1;
for j1=1:length(L11)
C1(L11(j1),:)=[1 0 0];
end
for j2=1:length(L12)
C1(L12(j2),:)=[0 0 1];
end
for j3=1:length(L21)
C2(L21(j3),:)=[1 0 0];
end
for j4=1:length(L22)
C2(L22(j4),:)=[0 0 1];
end
scatter(log2(Y_0(:,10)),log2(Y_0(:,4)),30,C2,'filled')
scatter(log2(Y_0(:,10)),log2(Y_0(:,4)),30,C1,'filled')