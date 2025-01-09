%% This MATLAB script generates the effect of OPN knockout on the tumor cells and IFNG
% Purpose of the study: To demonstrate OPN knockout improves the overall response to therapy 
y_0=[872.1;798.9;1090.2;1412.70;1253.4;1584.7;0;459.8;599.9;1474.4;1918.8;39.3;172.1;111.80];                                % Initial condition
alpha=[0.05 0.1 0.15 0.2 0.3 0.5];                                                                                           % Varrying immune accessibility
tspan=[0 70000]; 
anti_PD1=1;                                                                                                                  
OPN_Knockout=[0:10:500];                                                                                                     % OPN clearance rate
OPN_Knockout_max=500;
FC_Tum=[];
IFNG=[];
I=[];
for l=1:1:length(alpha)
P=ST_parameters_Fibro_Rich_OPN(alpha(l));                             
P(15)=1;                                                                                                                     % Setting spatial competition
[t_pre,x_pre]=ode15s(@(t,y)ST_mod_FIbro_Rich_OPN(t,y,0,P,0),[0 70000],y_0);
for k=1:1:length(OPN_Knockout)
[t_post,x_post]=ode15s(@(t,y)ST_mod_FIbro_Rich_OPN(t,y,anti_PD1,P,OPN_Knockout(k)),[0 70000],x_pre(length(x_pre),:));
x_post(x_post<1)=0;
% L2=find(x_post>10^4);
% x_post(L2,:)=[];
% t_post(L2)=[];
FC_Tum(k)=sum(x_post(length(x_post),1:4))/sum(x_pre(length(x_pre),1:4));
IFNG(k)=x_post(length(x_post),14)/(P(55)*P(31)/P(56));
end
subplot(1,2,1)
plot(OPN_Knockout/OPN_Knockout_max,FC_Tum,'LineWidth',2,'LineStyle','-')
xlabel('Concentration of OPN inhibitor','FontSize',18,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Fold Change (Tumor cells)// Pre vs Post ICI (FC_{X_C})','FontSize',18,'FontWeight','bold','FontName','Palatino Linotype')    % Plotting fold change of tumor cells vis-a-vis ICI wrt OPN clearance rate for different accessibility
hold on
subplot(1,2,2)
plot(OPN_Knockout/OPN_Knockout_max,IFNG,'LineWidth',2,'LineStyle','--')
xlabel('Concentration of OPN inhibitor','FontSize',18,'FontWeight','bold','FontName','Palatino Linotype')                            % Plotting IFNG concentration wrt OPN clearance rate across immune accessibilities
ylabel('IFN\gamma Levels (X_{IFN\gamma})','FontSize',18,'FontWeight','bold','FontName','Palatino Linotype')
I(l)=(1-tanh(alpha(l)*P(15)*x_pre(length(x_pre),9)))*100;
hold on
end
