%% This MATLAB code generates the main text plots for FIG. 1 of the solid tumor manuscript.
%% Figure 1 (A): Theoretical prediction of the relation between proliferation vs exhaution rate of the killer T cells for immune desert condition.
% Purpose of the study: To demonstarte for a given proliferation rate, there exists a threshold exhaustion rate for desertification

K_d=10;                  % Fixing the death rate of killer T cells
K_Tex=[0:1.5:25];        % Exhaustion rate
M2_1=[0:0.5:100];        % M_2 population
M2_2=[110:500:2500];
M2=[M2_1 M2_2];
K_prol=[];
for k=1:1:length(M2)
K_prol=K_Tex*M2(k)/(M2(k)+1)+K_d;                                            % Theoretical prediction for Proliferation-exhaustion relation for logistic growth.
plot(K_Tex,K_prol,'LineWidth',2);
set(gca,'FontSize',19,'FontWeight','bold','FontName','Palatino Linotype')
ylim([0 50])
xlabel('T cell exhaustion rate','FontSize',19.8,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('T cell proliferation rate','FontSize',19.8,'FontWeight','bold','FontName','Palatino Linotype')
hold on
box off
end

%% Figure 1(b-c): Computational validation of the theoretical result illustrated in figure 1
%% Figure 1(b): Simulation of killer vs exhausted state trajectories for different exhaustion rate with fixed proliferation rate
K_prol=20;                                                             % Fixing the proliferation rate
K_TKTex=[5 8 9.5 10.2 12 20];                                          % Simulating with different exhaustion rate
tspan=[0 70000];
y_0=[4092;1000;5000;0.02;0;1;300;100;10;0];                            % Initial condition          
P=ST_parameters_Immune_Desert;                                         % Nominal default parameters for immune desert
P(25)=10;                                                              % Death rate of Killer T cells    
P(26)=10;                                                              % Death rate of exhausted T cells
P(19)=15;                                                              % Proliferation of exhausted T cells 
P(17)=K_prol;                                                          % Proliferation of killer T cells
P(23)=0.5;                                                             % Low IL-2 based proliferation
figure
axes1 = axes('Position',[0.370277777777778 0.630838774485183 0.150295138888889 0.240973768588557]);
for k=1:1:length(K_TKTex)
P(21)=K_TKTex(k);
[t_pre,x_pre]=ode15s(@(t,y)ST_mod_Immune_Desert(t,y,0,0,0,P),tspan,y_0);
plot(log2(1+x_pre(:,6)/y_0(6)),log2(1+x_pre(:,4)/y_0(4)),'LineWidth',2)                                     % Plotting the killer vs the exhausted T cells trajectories for different exhaustion rates
xlabel('Exhausted T cell population','FontSize',19.2,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Killer T cell population','FontSize',19.2,'FontWeight','bold','FontName','Palatino Linotype')
hold on
set(axes1,'FontName','Palatino Linotype','FontSize',19,...
    'FontWeight','bold','LineWidth',1,'XTick',[0 2 4 6 8 10],'XTickLabel',...
    {'2^0','2^2','2^4','2^6','2^8','2^{10}'},'YTick',[0 2 6 10 14 18],'YTickLabel',...
    {'2^0','2^2','2^6','2^{10}','2^{14}','2^{18}'});
box off
end

%% Figure 1(c): Simulation of killer vs exhausted state trajectories for different proliferation rate with fixed exhaustion rate

K_TkTex=15;                                                             % Fixing the exhaution rate                                                            
P(21)=K_TkTex;
K_prol=[15 20 24 25 30 35];                                             % Simulating for different proliferation rate
figure
axes2 = axes('Position',...
    [0.370277777777778 0.630838774485183 0.150295138888889 0.240973768588557]);
for j=1:1:length(K_prol)
P(17)=K_prol(j);
[t_pre,x_pre]=ode15s(@(t,y)ST_mod_Immune_Desert(t,y,0,0,0,P),tspan,y_0);
plot(log2(1+x_pre(:,6)/y_0(6)),log2(1+x_pre(:,4)/y_0(4)),'LineWidth',2)                                    % Plotting the killer vs the exhausted T cells trajectories for different proliferation rates
xlabel('Exhausted T cell population','FontSize',19.2,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Killer T cell population','FontSize',19.2,'FontWeight','bold','FontName','Palatino Linotype')
hold on
set(axes2,'FontName','Palatino Linotype','FontSize',19,...
    'FontWeight','bold','LineWidth',1,'XTick',[0 2 4 6 8 10],'XTickLabel',...
    {'2^0','2^2','2^4','2^6','2^8','2^{10}'},'YTick',[0 2 6 10 14 18],'YTickLabel',...
    {'2^0','2^2','2^6','2^{10}','2^{14}','2^{18}'});
box off
end


%% Figure 1(d-f): Exploring the attractor landscape across the initial conditions.
% Purpose of the study: Uniquesness of an attractor governs the transition possiblity this is crucial for designing combination therapies for ICI non-responding attractors. 

y_0=5000*lhsdesign(300,12);                                                 % Generating a cohort of 300 different initial conditions within the carrying capacity
y_0(:,5)=0;                                                                 % Zero PD1- killer T cells in a pre-ICI scenario.    
y_1=5*lhsdesign(200,12);                                                    % Sampling the initial conditions near zero to ensure exhaustive sampling
y_1(:,5)=0;
Y_0=[y_0;y_1];
tspan=[0 70000];
P=ST_parameters_Immune_Rich;                                                % Nominal parameters for immune rich
P(17)=5;
P(20)=20;
P(21)=0.2;
P(22)=25;
P(24)=30;
X_CAF_T_DR=[];                                                              % For storing the final CAF and killer T cell components of the the attractors corresponding immune desert and rich phenotypes.   
ax1=axes();
for k=1:1:length(Y_0(:,1))
[t_pre,x_pre]=ode23s(@(t,y)ST_mod_Immune_Rich(t,y,0,0,P),tspan,Y_0(k,:));   
X_CAF_T_DR(k,:)=[x_pre(length(x_pre),8)/5000 x_pre(length(x_pre),4)/5000];
plot(x_pre(:,8)/5000,x_pre(:,4)/5000)                                       % Plotting CAF vs Killer T trajectory
hold on
scatter(X_CAF_T_DR(:,1),X_CAF_T_DR(:,2),300,[0 0 0],"x")
set(ax1,'FontName','Palatino Linotype','FontSize',19,'FontWeight','bold');
xlabel('CAF population','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Killer T population','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
hold on
box off
end

figure
P=ST_parameters_Immune_Rich;         
X_CAF_T_R=[];                                                               % For storing the final CAF and killer T cell components of the the attractors corresponding immune/fibro rich phenotypes. 
ax2=axes();
for k1=1:1:length(Y_0(:,1))
[t_pre,x_pre]=ode15s(@(t,y)ST_mod_Immune_Rich(t,y,0,0,P),[0 70000],Y_0(k1,:));
plot(x_pre(:,8)/(5000),x_pre(:,4)/(5000))
hold on
X_CAF_T_R(k1,:)=[x_pre(length(x_pre),8)/5000 x_pre(length(x_pre),4)/5000];
hold on
scatter(X_CAF_T_R(:,1),X_CAF_T_R(:,2),300,[0 0 0],"x")                                    
set(ax2,'FontName','Palatino Linotype','FontSize',19,'FontWeight','bold');
xlabel('CAF population','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')  
ylabel('Killer T population','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
hold on
box off
end

figure
P=ST_parameters_Immune_Rich;
P(17)=0.55;
P(20)=20;
P(21)=0.2;
P(22)=5;
P(24)=30;
X_CAF_T_D=[];                                                                 % For storing the final CAF and killer T cell components of the the attractors corresponding absolute immune desert phenotypes. 
ax3=axes();

for k2=1:1:length(Y_0(:,1))
[t_pre1,x_pre1]=ode15s(@(t,y)ST_mod_Immune_Rich(t,y,0,0,P),[0 70000],Y_0(k2,:));
if (x_pre1>-0.5 & x_pre1< 2*10^4)
plot(x_pre1(:,8)/(5000),x_pre1(:,4)/(5000))
hold on
X_CAF_T_D(k2,:)=[x_pre1(length(x_pre1),8)/5000 x_pre1(length(x_pre1),4)/5000];
hold on
scatter(X_CAF_T_D(:,1),X_CAF_T_D(:,2),300,[0 0 0],"x")
set(ax3,'FontName','Palatino Linotype','FontSize',19,'FontWeight','bold');
xlabel('CAF population','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Killer T population','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
hold on
box off
end
end


