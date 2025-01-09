%% This MATLAB code generates Figure 4 of the solid tumor manuscript. 
%% Figure 4(a-b) M2 knockout study
% Purpose of the study: To demonstrate that M2 knockout transition the TME
% trajectories from exhaustion-driven immune desert to Immune rich scenario. This is NOT true for low-proliferation induced immune desert

y_0=[4092;1000;5000;0.2;0;1;300;100;10;0];   % Initial condition
tspan=[0 70000];                             
anti_PD1=1;                                  % ICI drug input
IL2_level=[];                                 
u_IL2=1;                                     % External IL-2 injection

P=ST_parameters_Immune_Desert;               % Nominal parameter values for immune desert
P(25)=0.5;                                   % Death rate for killer T cells   
P(26)=3;                                     % Death rate for exhausted T cells
P(23)=0.5;                                   % Low IL-2 based proliferation
M2_Knock=[0 10^15];                          % M2 clearance rates
 
figure
subplot(1,2,1)
P(21)=25;                                                                                                     % Conversion rate from killer to exhausted T cells (Exhaustion-driven immune desert)
[t_pre_Exh,x_pre_Exh]=ode15s(@(t,y)ST_mod_Immune_Desert(t,y,0,0,M2_Knock(1),P),tspan,y_0);
[t_post_Knock_Exh,x_post_Knock_Exh]=ode15s(@(t,y)ST_mod_Immune_Desert(t,y,0,0,M2_Knock(2),P),tspan,y_0);      
plot(t_post_Knock_Exh,x_post_Knock_Exh(:,4)/y_0(4),'LineWidth',2,'LineStyle','-')                             % Plotting Killer T cells wrt time with M2 knockout
hold on
plot(t_pre_Exh,x_pre_Exh(:,4)/y_0(4),'LineWidth',2,'LineStyle','-')                                           % Plotting Killer T cells wrt time without M2 knockout
xlim([0 0.25])
xlabel('Time','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Killer T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')

P=ST_parameters_Immune_Desert;
P(17)=2;                                                                                                       % Low-proliferation induced immune desert                                                              
y_0(4)=10^-4;                                                                                                  % Setting the initial killer T cell population close to zero to ensure it falls in the domain of attraction
[t_pre_LP,x_pre_LP]=ode15s(@(t,y)ST_mod_Immune_Desert(t,y,0,0,M2_Knock(1),P),tspan,y_0);
[t_post_Knock_LP,x_post_Knock_LP]=ode15s(@(t,y)ST_mod_Immune_Desert(t,y,0,0,M2_Knock(2),P),tspan,y_0);
subplot(1,2,2)
plot(t_post_Knock_LP,x_post_Knock_LP(:,4)/y_0(4),'LineWidth',2,'LineStyle','-')                                % Plotting Killer T cells wrt time with M2 knockout
hold on   
plot(t_pre_LP,x_pre_LP(:,4)/y_0(4),'LineWidth',2,'LineStyle','-')                                              % Plotting Killer T cells wrt time with M2 knockout1
xlabel('Time','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Killer T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
xlim([0 0.3])

%% Figure 4(c-d) IL2 administer study
% Purpose of the study: To demonstrate that one time injection of a fibro-effected IL2 insertion can transition from an conditional desert to immune rich phenotype 

figure
P(24)=20;                                                                                      
F_IL2=[3:0.1:3.7];                                                                         % Different one-time IL-2 insertion amount
for k=1:1:length(F_IL2)
y_0(10)=F_IL2(k);                                                                          % Changing the IL-2 initial condition
[t_pre_LP,x_pre_LP]=ode15s(@(t,y)ST_mod_Immune_Desert(t,y,0,0,M2_Knock(1),P),tspan,y_0);
subplot(1,2,1)
plot(t_pre_LP,x_pre_LP(:,4)/y_0(4),'LineWidth',2,'LineStyle','-')                          % Plotting the time profile of the killer T cells for different initial IL-2
xlabel('Time','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Killer T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
hold on
xlim([0 3])
ylim([0 0.4*10^5])
end

P(28)=0.005;                                                                                % Ensuring Fibro-desert scenario
for k1=1:1:length(F_IL2) 
y_0(10)=F_IL2(k1);
[t_pre_LP,x_pre_LP]=ode15s(@(t,y)ST_mod_Immune_Desert(t,y,0,0,M2_Knock(1),P),tspan,y_0);
subplot(1,2,2)
plot(t_pre_LP,x_pre_LP(:,4)/y_0(4),'LineWidth',2,'LineStyle','-')
xlabel('Time','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')             % Plotting the time profile of the killer T cells for different initial IL-2 in fibro-desert situation
ylabel('Killer T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
hold on
xlim([0 3])
ylim([0 0.4*10^5])
end

%% Figure 4(e-f) OPN knockout study
% Purpose of the study: To demonstrate that one time injection of a fibro-effected IL2 insertion can transition from an conditional desert to immune rich phenotype 

figure
y_0=[4000;1000;10;10;5000;50;0;0;300;100;10;0;10;10];                                            % Initial condition
OPN_Knockout1=[500:500:2500];                                                                    % Different OPN clearance rates
alpha=0.8;                                                                                       % Very low immune accessibility
P=ST_parameters_Fibro_Rich_OPN(alpha);
for k=1:1:length(OPN_Knockout1)
[t_pre,x_pre]=ode15s(@(t,y)ST_mod_FIbro_Rich_OPN(t,y,0,P,OPN_Knockout1(k)),[0 70000],y_0);
subplot(1,2,1)
plot(t_pre,x_pre(:,9)/5000,'LineWidth',2)                                                                         % Plotting the time profiles of CAF for different OPN clearnace rate
xlim([0 2])
xlabel('Time','FontSize',19,'FontWeight','bold','FontName','Palatino Linotype')                          
ylabel('CAF population\\(Pre-ICI)','FontSize',19,'FontWeight','bold','FontName','Palatino Linotype')                  
hold on
end
y_0(9)=2000;
subplot(1,2,2)
[t_pre,x_pre]=ode15s(@(t,y)ST_mod_FIbro_Rich_OPN(t,y,0,P,OPN_Knockout1(1)),[0 70000],y_0);
plot(x_pre(1:200,6)/5000,x_pre(1:200,9)/5000,'LineWidth',2,'LineStyle','--')
xlabel('Killer T cell population//(Pre-ICI)','FontSize',19,'FontWeight','bold','FontName','Palatino Linotype')      % Plotting CAF vs. killer T with and without OPN.   
ylabel('CAF population\\(Pre-ICI)','FontSize',19,'FontWeight','bold','FontName','Palatino Linotype')
hold on
[t_pre,x_pre]=ode15s(@(t,y)ST_mod_FIbro_Rich_OPN(t,y,0,P,OPN_Knockout1(end)),[0 70000],y_0);
plot(x_pre(1:350,6)/5000,x_pre(1:350,9)/5000,'LineWidth',2,'LineStyle','-')


% K_TProl=[0.5:0.5:30];
% K_CAFProl=[0.05:(9-0.05)/60:9-(9-0.05)/60];
% CAF_ST=[];
% TK_ST=[];
% P=ST_parameters_Immune_Desert;
% %P(25)=10;
% P(24)=0.02;
% P(29)=1.5;
% P(30)=1.5;
% P(31)=10;
% t_span=[0 70000];
% y_0=[4092;1000;5000;0.2;0;1;30;10;10;0];
% for k=1:1:length(K_TProl)
%         P(17)=K_TProl(k);
%     for j=1:1:length(K_CAFProl)
% P(28)=K_CAFProl(j);
% [t_pre,x_pre]=ode15s(@(t,y)ST_mod_Immune_Desert(t,y,0,0,0,P),[0 70000],y_0);
% CAF_ST(k,j)=x_pre(length(x_pre),7)/5000;
% TK_ST(k,j)=x_pre(length(x_pre),4)/5000;
% CAF_ST(CAF_ST<0.0005)=0;
% TK_ST(TK_ST<0.0005)=0;
% end
% end
% CAF_Dum=zeros(length(CAF_ST));
% P=imfuse(CAF_ST,CAF_Dum);
% P_CAF=P;
% P_CAF(:,:,1)=P(:,:,2);
% P_CAF(:,:,2)=P(:,:,1);
% 
% TK_Dum=zeros(length(TK_ST));
% P1=imfuse(TK_ST,TK_Dum);
% P_TK=P1;
% P_TK(:,:,3)=P1(:,:,2);
% P_TK(:,:,2)=P1(:,:,3);
% 
% P_Tot=(P_CAF+P_TK);
% imshow(P_Tot)