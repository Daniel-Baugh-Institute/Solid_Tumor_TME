%% This MATLAB script generates the figures corresponding to the fibro-desert phenotype in the supplementary figure
% Purpose of the study: A higher exhausted to tumor cell interaction increases PDL1-/PDL1+ ratio but the post-ICI interaction scenario remains unaffected.


%% Figure (a): Pre-ICI scenario
K=1200;                                        
y_0=[1000;1000;1000;5000;50;0;10];                                           % Initial condition
P=ST_parameters_Fibro_desert;                                                % Nominal parameter values for fibro-desert
P1=P; 
t_end=50000;
x_ST=[];
Counter=1:10:K;
T1=[20 400 800 1600 3200 6400 12800 20000];                                 % Different Exhauted T cell-driven growth rate of tumor cells
ax=axes();
for k=1:length(T1)                                                                                                   % Loading exhausted T cell-driven growth
    P1(9)=T1(k);
    [t_pre,x_pre]=ode15s(@(t,y)ST_mod_FIbro_Desert(t,y,0,P1),[0 t_end],y_0);
x_ST(k,:)=x_pre(length(x_pre(:,1)),:);
% Conditioning the solution space
x_pre(x_pre(:,1)<0)=0;
x_pre(x_pre(:,2)<0)=0;
plot((x_pre(:,1)+x_pre(:,2))/30000,x_pre(:,3)/30000,'LineWidth',2,'LineStyle','--')                                  % Plotting total PDL1- fraction vs. PDL1+ fraction
set(ax,'FontName','Palatino Linotype','FontSize',19,'FontWeight','bold','LineWidth',1);
xlabel('Non-PDL1 tumor cell population','FontSize',21,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('PDL1^+ tumor cell population','FontSize',21,'FontWeight','bold','FontName','Palatino Linotype')
box off
hold on 
end


%% Figure (b): Post-ICI scenario
figure
K=151;
y_0=[1000;1000;1000;5000;50;0;10];                                          % Initial condition
P=ST_parameters_Fibro_desert;                                               % Nominal parameters for fibro-desert   
P1=P;
t_end=50000;
x_ST=[];
x_ST1=[];
Counter=1:50:K;                                           
anti_PD1=2;                                                                 % Seting anti-PD1 level
T1=[20 400 800 1600 3200 6400 12800 25000];                                 % Tumor cell growth by exhausted T cells
for k=1:length(T1)
P1(9)=T1(k);                                                                        % Modulating the exhausted T cell-driven proliferation of tumor cells                                                            
[t_pre,x_pre]=ode15s(@(t,y)ST_mod_FIbro_Desert(t,y,0,P1),[0 t_end],y_0);           
x_ST(k,:)=x_pre(length(x_pre(:,1)),:);
x_pre(x_pre(:,1)<0)=0;                                                          
x_pre(x_pre(:,2)<0)=0;
[t_post,x_post]=ode15s(@(t,y)ST_mod_FIbro_Desert(t,y,anti_PD1,P1),[0 t_end],y_0); 
x_ST1(k,:)=x_post(length(x_post(:,1)),:);                                   
x_post(x_post(:,1)<0)=0;
x_post(x_post(:,2)<0)=0;
subplot(1,3,1)
plot(x_post(:,6)/5000,(x_post(:,1)+x_post(:,2)+x_post(:,3))/30000,'LineWidth',2,'LineStyle','-')   % Plotting PD1- killer T cells vs total post-ICI tumor cell population
hold on
box off
subplot(1,3,2)
plot(x_pre(:,5)/5000,(x_pre(:,1)+x_pre(:,2)+x_pre(:,3))/30000,'LineWidth',2,'LineStyle','--')      % Plotting PD1+ killer T cells vs total pre-ICI tumor cell population
hold on
box off
subplot(1,3,3)
plot(t_post,(x_post(:,1)+x_post(:,2)+x_post(:,3))/30000,'LineWidth',2,'LineStyle','-')             % Plotting total post-ICI tumor cell population wrt time
xlim([0 1.75])
xlabel('Time','FontSize',19,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Total tumor cells','FontSize',19,'FontWeight','bold','FontName','Palatino Linotype')
hold on
box off
end
