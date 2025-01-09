%% This MATLAB script generates the plots for immune-rich solid tumor subtype

%% Initialisation
anti_PD1=1;
tspan=[0 70000];

%% Tumor cells
C_NPDL10=9990;
C_PDL10=8500;
C_0=[C_NPDL10;C_PDL10];

%% Resources
R_0=12;

%% T cells
T_KPD10=2;
T_KNPD10=0;
TEX0=65;
T_0=[T_KPD10;T_KNPD10;TEX0];

%% CAF 
FWT_0=100;
CAF_0=4500;

%% Macrophages
M_10=84;
M_20=48;
M_0=[M_10;M_20];

%% Interlukin-8
IL8_0=25;

%%  IL2
IL2_0=0;

y_0=[C_0;R_0;T_0;FWT_0;CAF_0;M_0;IL8_0; IL2_0];
u_Tk=[0:5:5000];                                    % External insertion of killer T cells
P=ST_parameters_Immune_Rich;
CAF_ST=[];
TK_ST=[];
%% Fig CAF vs Killer T cells
% Purpose of the study: Killer T cells slightly increase the CAF population in the neighborhood of an immune rich attractor
for k=1:1:length(u_Tk)
    [t_pre, x_pre]=ode15s(@(t,y)ST_mod_Immune_Rich(t,y,0,u_Tk(k),P),tspan,y_0);
    CAF_ST(k)=x_pre(length(x_pre(:,1)),8);
    TK_ST(k)=x_pre(length(x_pre(:,1)),4);
    plot(x_pre(:,4)/5000,x_pre(:,8)/5000,'LineWidth',2)                                                         % Plotting CAF vs the Killer T cells                        
    xlabel('PD1^+ Cytotoxic T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
    ylabel('CAF','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
    hold on
end

%% Fig B: Post-ICI scenario of CAF, Killer T, and Total tumor cells
% Purpose of the study: CAF population reduces post-ICI in immune-rich attractor
figure
P=ST_parameters_Immune_Rich;
P(3)=P(1)/10;                                                                                                     % Low PDL1+ natural proliferation
Q=P;
CAF_FC=[];
CAF_prol=[0.005:0.16:0.325]
for m=1:1:length(CAF_prol)
    Q(32)=P(32)*CAF_prol(m);
[t_pre, x_pre]=ode15s(@(t,y)ST_mod_Immune_Rich(t,y,0,0,Q),tspan,y_0);
    x_pre(x_pre<0)=0;
    [t_post, x_post]=ode15s(@(t,y)ST_mod_Immune_Rich(t,y,anti_PD1,0,Q),tspan,y_0+0.005);
    CAF_ST(m)=x_pre(length(x_pre(:,1)),8)/x_post(length(x_post(:,1)),8);
    plot3((x_pre(:,8))/5000,(x_pre(:,4))/5000,(x_pre(:,1)+x_pre(:,2))/20000,'LineWidth',2,'LineStyle','--')        % Plotting pre-ICI tumor vs killer T vs. CAF cells
    hold on
    plot3((x_post(:,8))/5000,(x_post(:,5))/5000,(x_post(:,1)+x_post(:,2))/20000,'LineWidth',2,'LineStyle','-')     % Plotting post-ICI tumor vs killer T vs. CAF cells
    xlabel('CAF','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
    ylabel('Killer T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
    zlabel('Tumor cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
hold on
end

%% Fig C: Pre and Post ICI attractor landscape for immune rich
% Purpose of the study: To demonstarte the pre-ICI and post-ICI attractor differs in increased killer T cells and decreased CAF population
y_0=5000*lhsdesign(1000,12);                                                % Generating 1000 initial conditions
y_0(:,5)=0;                                                                 % Ensuring zero PD1- killer T cells before ICI
Y_0=[y_0];
tspan=[0 70000];                                           
P=ST_parameters_Immune_Rich;
P(3)=P(1)/10;
X_CAF_T_DR=[];
X_CAF_T_DR_post=[];
ax1=axes();
for k=1:1:length(Y_0(:,1))
[t_pre,x_pre]=ode23s(@(t,y)ST_mod_Immune_Rich(t,y,0,0,P),tspan,Y_0(k,:));
X_CAF_T_DR(k,:)=[x_pre(length(x_pre),8)/5000 x_pre(length(x_pre),4)/5000];
plot(x_pre(:,8)/5000,x_pre(:,4)/5000)                                                             % Plotting Killer T vs. CAF cells pre-ICI
hold on
scatter(X_CAF_T_DR(:,1),X_CAF_T_DR(:,2),300,[0 0 0],"x")                                          % Plotting the pre-ICI attractor points
set(ax1,'FontName','Palatino Linotype','FontSize',19,'FontWeight','normal');
xlabel('CAF population','FontSize',20,'FontWeight','normal','FontName','Palatino Linotype')
ylabel('Killer T population','FontSize',20,'FontWeight','normal','FontName','Palatino Linotype')
hold on
box off
end

figure
X_CAF_T_DR_post=[];
ax2=axes();
for k=1:1:length(Y_0(:,1))
[t_post,x_post]=ode23s(@(t,y)ST_mod_Immune_Rich(t,y,2,0,P),tspan,Y_0(k,:));                    
X_CAF_T_DR_post(k,:)=[x_post(length(x_post),8)/5000 x_post(length(x_post),5)/5000];              % Plotting Killer T vs. CAF cells post-ICI
L1=find(x_post(:,5)>5005);                                                                             
x_post(L1,:)=[];t_post(L1)=[];
plot(x_post(:,8)/5000,x_post(:,5)/5000)
hold on
scatter(X_CAF_T_DR_post(:,1),X_CAF_T_DR_post(:,2),300,[0 0 0],"x")                               % Plotting the pre-ICI attractor points

set(ax2,'FontName','Palatino Linotype','FontSize',19,'FontWeight','normal');
xlabel('CAF population','FontSize',20,'FontWeight','normal','FontName','Palatino Linotype')
ylabel('Killer T population','FontSize',20,'FontWeight','normal','FontName','Palatino Linotype')
xlim([0 1])
ylim([0 1])
hold on
box off
end

% figure
% %CAF_C=[10:120:6*10^4];
% CAF_C=[240 720 1200 1680 2160];
% P=ST_parameters_Immune_Rich;
% P(3)=P(1)/10;
% Q=P;
% for l=1:1:length(CAF_C)
%     Q(10)=CAF_C(l)*P(11);
%     [t_pre, x_pre]=ode15s(@(t,y)ST_mod_Immune_Rich(t,y,0,0,Q),tspan,y_0);
%     x_pre(x_pre<0)=0;
%     x_pre(x_pre>10^4)=10^4;
%     subplot(1,2,2)
%     plot((x_pre(:,8))/5000,x_pre(:,1)/20000,'LineWidth',2)
%     xlabel('CAF','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
%     ylabel('PDL1^- Tumor cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
%     hold on
% end

% figure
% P=ST_parameters_Immune_Rich;
% P(3)=P(1)/10;
% Q=P;
% CAF_C=[0:60:2400];
% CAF_Prol1=[0.0005:0.005:0.5];
% CAF_Prol2=[0.5:15:60];
% CAF_Prol=[CAF_Prol1 CAF_Prol2];
% Tum_Tot=[];
% for n=1:1:length(CAF_Prol)
% Q(32)=CAF_Prol(n);
%     for l=1:1:length(CAF_C)
%     Q(10)=CAF_C(l)*P(11);
%     [t_pre, x_pre]=ode15s(@(t,y)ST_mod_Immune_Rich(t,y,0,0,Q),tspan,y_0);
%     x_pre(x_pre<0)=0;
%     x_pre(x_pre>10^4)=10^4;
%     Tum_Tot(l)=(x_pre(length(x_pre),1)+x_pre(length(x_pre),2))/20000;
%     end
%     plot(CAF_C,Tum_Tot,'LineWidth',2)
%     xlabel('\zeta','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
%     ylabel('Total tumor cell population','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
%     hold on
%     box off
% end