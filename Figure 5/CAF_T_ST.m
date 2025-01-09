%% Fig C: Pre and Post ICI attractor landscape for immune rich
% Purpose of the study: To demonstarte the pre-ICI and post-ICI attractor differs in increased killer T cells and decreased CAF population
y_0=5000*lhsdesign(1,12);                                                % Generating 1000 initial conditions
y_0(:,5)=0;                                                                 % Ensuring zero PD1- killer T cells before ICI
Y_0=[y_0];
tspan=[0 70000];                                           
P=ST_parameters_Immune_Rich;
P(17)=0.5;
P(18)=0.5;
P(21)=100;
P(23)=150;
P(1)=2*P(1);
P(3)=P(1)/10;
P(51)=0.002;
anti_PD1=2;
Y_0=5000*lhsdesign(1000,12);
Y_0(:,5)=0;
X_ST=[];
ax2=axes();
P(11)=5000;
for k=1:1:1000
[t_pre,x_pre]=ode23s(@(t,y)ST_mod_Immune_Rich(t,y,anti_PD1,0,P),tspan,Y_0(k,:)');
X_ST(k,:)=x_pre(length(x_pre),:);
plot(t_pre,log2(x_pre(:,1)+x_pre(:,2)),'LineWidth',2)
set(ax2,'FontName','Palatino Linotype','FontSize',19,'FontWeight','normal');
xlabel('Time','FontSize',20,'FontWeight','normal','FontName','Palatino Linotype')
ylabel('Post-ICI Tumor cells','FontSize',20,'FontWeight','normal','FontName','Palatino Linotype')
xlim([0 25])
hold on
box off
end
C=zeros(1000,3);
L1=find(log2(X_ST(:,1)+X_ST(:,2))<10);
L2=find(log2(X_ST(:,1)+X_ST(:,2))>10);
for k1=1:1:length(L1)
    C(L1(k1),:)=[1 0 0];
end
for k2=1:1:length(L2)
    C(L2(k2),:)=[0 0 1];
end
figure
% ax3=axes();
% scatter(log2(X_ST(:,5)),log2(X_ST(:,8)),50,C,'filled')
% set(ax3,'FontName','Palatino Linotype','FontSize',19,'FontWeight','normal');
% xlabel('Killer T cell','FontSize',20,'FontWeight','normal','FontName','Palatino Linotype')
% ylabel('CAF population','FontSize',20,'FontWeight','normal','FontName','Palatino Linotype')

TF=[X_ST(:,5) X_ST(:,8)];
TF=round(TF*100)/100;
U1=unique(TF(:,1));
U2=unique(TF(:,2));
S=[];
for j=1:1:length(U1)+length(U2)
    clust=kmeans(TF,j);
    s= silhouette(TF,clust,'cityblock');
S(j)=mean(s);
end
Lf=find(S==max(S));
[Clust,C]=kmeans(TF,Lf)
subplot(1,2,1)
bar(2,log2(C(1,:)))
subplot(1,2,2)
bar(2,log2(C(2,:)))