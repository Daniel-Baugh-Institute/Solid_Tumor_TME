%% This MATLAB code generates the necessary figures for Figure 1 of the manuscript

%% Figure 1B
%% Purpose of the study: Characterizing the attractor landscape for the core solid tumor network
% Load the parameter file P_para_ST_FIG

alpha=[0.01 0.15 0.25]; % Different immune accessibility
y_0=10*lhsdesign(12,1); % Fix the initial condition
y_0(7)=0;
y_0=[6.7996 0.3811 6.5338 2.1701 9.9400 4.8818 0 8.2925 3.2460 1.4985 8.7039 3.7097];
 %y_0 =[10;5;10;5;1;2;0;2;3;4;4;2];
y_0=[6.7996 3.3811 6.5338 5.1701 9.9400 4.8818 0 8.2925 3.2460 1.4985 8.7039 3.7097];

tspan=[0 70000]; 
X_ST=[]; % Contains the attractors for different parameter choices.
l=1;
S=[];
for j=1:1:length(P_para(1,:))
for k=1:1:length(alpha)
    P=ST_parameters(alpha(k));              % Loads the nominal parameter values for the model
        P(1:4)=P(1:4)/3;                    % Tuning the natural proliferation rates of the tumor cells
        P(9)=P_para(1,j);                   % Exhausted T cell driven growth of the tumor cells 
        P(10)=P_para(2,j);                  % CAF-driven growth of immune-accessible tumor cells 
        P(11)=1.5*P(10);                    % CAF-driven growth of immune-inaccessible tumor cells
        P(22)=2.5*P_para(3,j);              % Proliferation rate of PD1- T cells
        P(24)=P_para(11,j);                 % Exhaustion rate of killer T cells
        P(25)=P_para(4,j);                  % MHC sensing of tumor cells for immune activation
        P(26)=P_para(5,j);                  % CAF-driven inhibition of T cells through regulatory T cells.
        P(27)=P_para(6,j);                  % Death rate for killer PD1+ T cells 
        P(31)=P_para(7,j);                  % Tumor cells-induced growth of CAF
        P(32)=P_para(8,j);                  % Proximity factor for CAF and resistant tumor cells  
        P(33)=P_para(9,j);                  % M2-driven proliferation rate for CAFs 
        P(34)=P_para(10,j);                 % Death rate for CAF                         
[t_pre,x_pre]=ode23s(@(t,y)ST_mod(t,y,0,P),tspan,y_0);

% Eliminating infeasible solutions
if(x_pre<-0.5) 
    x_pre=[];
elseif (x_pre>10^5) 
    x_pre=[];
else

X_ST(l,:)=[log10(1+(x_pre(length(x_pre),1)+x_pre(length(x_pre),3))/(x_pre(length(x_pre),2)+x_pre(length(x_pre),4))) x_pre(length(x_pre),6)/5000 x_pre(length(x_pre),9)/5000];   % Storing the attractor values for PDL1-/PDL1+ tumor cell proportion, killer T, and CAF,
Iacc=tanh(alpha(k)*x_pre(length(x_pre),9)*P(14));    % Computing the immune accessibility index for different parameter scenario                                                     
Acc=1-Iacc;
S(l,:)=[Iacc 0 Acc];                                 % Color coding based on immune accessiblity
S(S<0)=0;
S(S>1)=1;
if (X_ST(l,1)<=0.05 && X_ST(l,2)<=0.05)
X_ST(l,:)=[];
x_pre=[];
S(l,:)=[];
else
l=l+1;
plot3(log10(1+(x_pre(:,1)+x_pre(:,3))./(x_pre(:,2)+x_pre(:,4))),x_pre(:,6)/5000,x_pre(:,9)/5000,'LineWidth',0.5,'Color',[0 0 0]); % Plotting different state trajectories for attractors.
%xlim([0 5])
ylim([-0.005 1])
zlim([-0.005 1])
hold on
%scatter3(X_ST(l,1),X_ST(l,2),X_ST(l,3),90,S(l,:),"x")
%hold on
end
end
end
end
hold on
scatter3(X_ST(:,1),X_ST(:,2),X_ST(:,3),90,S,"x")                            % Plotting the attractors on the state trajectories
P_para_Fibro_des=5*lhsdesign(4,10);                                         % Parameters for the fibro-desert scenario
alpha1=0.3;
X_ST_des=[];
S1=[];
p=1;
P=ST_parameters(alpha1);
for l1=1:1:length(P_para_Fibro_des(1,:))
        P(30)=P_para_Fibro_des(1,l1);
        P(31)=P_para_Fibro_des(2,l1);
        P(32)=P_para_Fibro_des(3,l1);
        P(33)=P_para_Fibro_des(4,l1);
        P(34)=max(P_para(:,25));
[t_pre_des,x_pre_des]=ode23s(@(t,y)ST_mod(t,y,0,P),tspan,y_0);
if(x_pre_des<-0.5)
    x_pre_des=[];
else
   plot3(log10(1+(x_pre_des(:,1)+x_pre_des(:,3))./(x_pre_des(:,2)+x_pre_des(:,4))),x_pre_des(:,6)/5000,x_pre_des(:,9)/5000,'LineWidth',0.5,'Color',[0 0 0]);
   xlim([0 3])
   ylim([-0.005 1])
   zlim([-0.005 1])
   hold on
Iacc=tanh(alpha1*x_pre_des(length(x_pre_des),9)*P(14));
 Acc=1-Iacc;
 S1(p,:)=[Iacc 0 Acc];
 S1(S1<0)=0;
 S1(S1>1)=1;
scatter3(log10(1+(x_pre_des(length(x_pre_des),1)+x_pre_des(length(x_pre_des),3))/(x_pre_des(length(x_pre_des),2)+x_pre_des(length(x_pre_des),4))),x_pre_des(length(x_pre_des),6)/5000,x_pre_des(length(x_pre_des),9)/5000,90,S1(p,:),"x")
hold on
X_ST_des(p,:)=[log10(1+(x_pre_des(length(x_pre_des),1)+x_pre_des(length(x_pre_des),3))/(x_pre_des(length(x_pre_des),2)+x_pre_des(length(x_pre_des),4))) x_pre_des(length(x_pre_des),6)/5000 x_pre_des(length(x_pre_des),9)/5000];
p=p+1;
end
end

