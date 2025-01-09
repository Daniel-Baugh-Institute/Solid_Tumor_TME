%% This MATLAB script generates the necessary plots for Fibro-rich scenario

%% Fig Il8 vs killer T cells
% Purpose of the study: To demonstarte IL8 increase post-ICI in fibro-rich conditions
%y_0=[20;1000;1000;1000;1000;30;0;10;10;10;0;5*10^4];
y_0=500*lhsdesign(12,1);                                                    % using any initial condition at random
y_0(7)=0;                                                                   % ensuring zer0o PD1- killer T cells during pre-ICI
x_ST1=[];                                                          
anti_PD1=1;                                                                 % fixing anti-PD1 level
tspan=[0 70000];
alpha=[0.005 0.5];                                                          % Immune accessible and inaccessible
I=zeros(length(alpha),1);
Path_TIL8_Fib_Rich=[];
for j=1:1:length(alpha)
P=ST_parameters_Fibro_Rich(alpha(j));
Y_T=P(21);
Y_IL8=(P(45)*P(36)+P(46)*P(30)+2*P(47)*P(7))/P(48);                                                % Calculating the carrying capacity of IL-8
[t_post_Fib_Rich,x_post_Fib_Rich]=ode15s(@(t,y)ST_mod_FIbro_Rich(t,y,anti_PD1,P),tspan,y_0);
[t_pre_Fib_Rich,x_pre_Fib_Rich]=ode15s(@(t,y)ST_mod_FIbro_Rich(t,y,0,P),tspan,y_0);
plot(x_pre_Fib_Rich(:,6)/Y_T,x_pre_Fib_Rich(:,12)/Y_IL8,'LineStyle','--','LineWidth',2);                 % Plotting pre-ICI IL8 concentration vs. killer T cell population
hold on
plot(x_post_Fib_Rich(:,7)/Y_T,x_post_Fib_Rich(:,12)/Y_IL8,'LineStyle','-','LineWidth',2);                % Plotting post-ICI IL8 concentration vs. killer T cell population
hold on
end


%% Fig Immune accessible vs immune inaccessible tumor cells
% Purpose of the study: At lower immune-inaccessibility the total tumor cell becomes unreachable from anti-PD1

% Finding the pre-ICI fibro-rich attractors for different immune accessiblity
alpha=[0.0005 0.05 0.1 0.3 0.8];
Y_0=5000*lhsdesign(300,12);
Y_0(:,7)=0;
figure
ax2=axes();
L1=[];
for k=1:1:length(Y_0(:,1))
for m=1:1:length(alpha)
    P=ST_parameters_Fibro_Rich(alpha(m));
    P(14)=0.001;                                                                                                                 % Low spatial competition across the tumor cells
    P(3)=P(1);
    [t_pre, x_pre]=ode15s(@(t,y)ST_mod_FIbro_Rich(t,y,0,P),[0 70000],Y_0(k,:)');
    L1(k,:,m)=x_pre(length(x_pre),:);
    plot((x_pre(:,1)+x_pre(:,2))/20000,(x_pre(:,3)+x_pre(:,4))/20000)                                                            % Plotting pre-ICI accessible vs inaccessible tumor cells                             
    hold on
    scatter((L1(k,1,m)+L1(k,2,m))/20000,(L1(k,3,m)+L1(k,4,m))/20000, 300,[0 0 0],"x")                                            % Superimposing the attractors
    hold on
    set(ax2,'FontName','Palatino Linotype','FontSize',19,'FontWeight','normal','LineWidth',1);
    xlabel('Immune-accessible tumor cell population','FontSize',21,'FontWeight','bold','FontName','Palatino Linotype')
    ylabel('CAF-protected tumor cell population','FontSize',21,'FontWeight','bold','FontName','Palatino Linotype')
    xlim([0 1])
    ylim([0 1])
    box off
    hold on
end
end

% The post ICI scenario
figure
ax3=axes();
for m1=1:1:length(alpha)
    P=ST_parameters_Fibro_Rich(alpha(m1));
    P(14)=0.001;
    P(3)=P(1);
    [t_pre, x_pre]=ode15s(@(t,y)ST_mod_FIbro_Rich(t,y,0,P),[0 70000],Y_0(1,:)');
%    [t_post, x_post]=ode15s(@(t,y)ST_mod_FIbro_Rich(t,y,1,P),[0 70000],L1(1,:,m1));
    [t_post, x_post]=ode15s(@(t,y)ST_mod_FIbro_Rich(t,y,1,P),[0 70000],x_pre(length(x_pre),:));
    plot((x_post(:,1)+x_post(:,2))/20000,(x_post(:,3)+x_post(:,4))/20000,'LineWidth',2)                                          % Plotting pre-ICI accessible vs inaccessible tumor cells 
     set(ax3,'FontName','Palatino Linotype','FontSize',19,'FontWeight','bold','LineWidth',1);
    xlabel('Immune-accessible tumor cell population','FontSize',21,'FontWeight','bold','FontName','Palatino Linotype')
    ylabel('CAF-protected tumor cell population','FontSize',21,'FontWeight','bold','FontName','Palatino Linotype')
    xlim([0 1])
    ylim([0 1])
    box off
    hold on
end