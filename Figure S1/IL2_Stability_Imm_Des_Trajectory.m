%% This MATLAB script generates the figures related to IL-2 modulation
% Purpose of the study: IL2 input fixed rate modulates the stability and one time IL-2 facilitates transition to immune rich attractor
y_0=[4000;1000;5000;0.2;0;0;300;100;10;0];                                  % Initial condition
tspan=[0 70000];
anti_PD1=1;
spec_A_pre=[];                                                   
spec_A_post=[];
A=[];
B=[];
IL2_level=[];
u_IL2=0.8;
P=ST_parameters_Immune_Desert;
P(23)=8;                                                                    % Tuning IL2-driven killer T cell proliferation rate
F_IL2=[0:0.001:3];                                                          % Setting different IL-2 initial spikes
Eig_pre=[];
Eig_post=[];
for k=1:1:length(F_IL2)
Eig_pre(k)=P(17)*(1+P(23)*(F_IL2(k)/P(42))/(F_IL2(k)/P(42)+1))-P(21)-P(25);    % Finding the spectral radius associated with the immune-desert attractor Il2 followed by ICI scenario
Eig_post(k)=P(18)*(1+P(23)*(F_IL2(k)/P(42))/(F_IL2(k)/P(42)+1))-P(21)-P(25);    % Finding the spectral radius associated with the immune-desert attractor ICI followed by IL2 scenario
end
ax1=axes();
plot(F_IL2,Eig_pre,'LineWidth',2)                                                                              % Plotting the spectral radius IL2 followed by ICI
hold on
plot(F_IL2,Eig_post,'Linewidth',2)                                                                             % Plotting the spectral radius ICI + IL2
xlabel('External IL-2 level','FontSize',21,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Stability margin (max(\sigma_A))','FontSize',21,'FontWeight','bold','FontName','Palatino Linotype')
set(ax1,'FontName','Palatino Linotype','FontSize',19,'FontWeight','bold','LineWidth',1);
box off


figure

ax2 = axes();
x_post(x_post<0)=0;
x_pre(x_pre<0)=0;
x_post(1,4)=2*x_pre(1,4);
plot(log((x_post(:,4)+x_post(:,5))/10000),log((x_post(:,1)+x_post(:,2))/20000),'LineWidth',2,'Color',[0 0 1])               % Plotting the Killer T vs Tumor in the IL-2 followed by ICI
xlabel('Killer T cell population (x_{T_k}(u_{IL2}^s))','FontSize',18,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Total tumor cell population (x_{Tum}(u_{IL2}^s))','FontSize',18,'FontWeight','bold','FontName','Palatino Linotype') 
hold on
plot(log((x_pre(:,4)+x_pre(:,5))/5000),log((x_pre(:,1)+x_pre(:,2))/20000),'LineWidth',2,'Color',[1 0 0])                    % Plotting the Killer T vs Tumor in the IL-2 + ICI
xlim([-20 0])
set(ax2,'FontName','Palatino Linotype','FontSize',19,'FontWeight','bold','LineWidth',1);
box off

figure
ax3 = axes();
plot((x_pre(:,4)+x_pre(:,5))/5000,(x_pre(:,1)+x_pre(:,2))/20000,'LineWidth',2,'Color',[1 0 0])
xlim([0 3.95*10^-5])
xlabel('Killer T cell population (x_{T_k}(u_{IL2}^s))','FontSize',18,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Total tumor cell population (x_{Tum}(u_{IL2}^s))','FontSize',18,'FontWeight','bold','FontName','Palatino Linotype')
set(ax3,'FontName','Palatino Linotype','FontSize',19,'FontWeight','bold','LineWidth',1);
box off
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';