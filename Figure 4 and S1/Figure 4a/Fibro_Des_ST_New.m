%% This MATLAB script generates the plots for the fibro-desert subtype figure 3 of the main manuscript
% Purpose of the study: (a) Existence of a scenario of high immune activity as defined in the manuscript
%(b) Resource supply rate governs the existence of the high immune activity area

Res=[10:1000:10^7];                                                          % Different resource rates                
y_0=[1000;1000;1000;5000;50;0;10];                                           % Initial condition
t_end=70000;                                                                 % Time span
X_ST=[];                                                                     % Stores pre-ICI attractor states
X_ST_post=[];                                                                % Stores post-ICI attractor states   
for k=1:length(Res)
    L_pre=[];K_pre=[]; L_post=[];K_post=[];
    P=ST_parameters_Fibro_desert;                                            % Nominal parameter values for fibro desert
    P(10)=30;                                                                % Specifying a higher cytotoxicity
    P(15)=0;                                                                 % No spatial competition
    P(8)=Res(k);                                                             % Modulating the maximum resource withholding capacity
    [t_pre,x_pre]=ode15s(@(t,y)ST_mod_FIbro_Desert(t,y,0,P),[0 t_end],y_0);
   % Conditioning the pre-ICI solution space
    [L_pre,K_pre]=find(x_pre(:,1:3)>10^4);
    x_pre(unique(L_pre),:)=[];
    t_pre(unique(L_pre))=[];
    X_ST(k,:)=x_pre(length(x_pre),:);                                                  % Loading the pre-ICI attractor set

   [t_post,x_post]=ode15s(@(t,y)ST_mod_FIbro_Desert(t,y,2,P),[0 t_end],X_ST(k,:)');
   % Conditioning the post-ICI solution space
   [L_post,K_post]=find(x_post(:,1:3)>10^4);
    x_post(unique(L_post),:)=[];
    t_post(unique(L_post))=[];
    X_ST_post(k,:)=x_post(length(x_post),:);                                            % Loading the post-ICI attractor set                                        
end


ax1=axes();
C=zeros(length(Res),3);
C(:,1)=0;
C(:,2)=(X_ST(:,1)+X_ST(:,2))/20000;
C(:,3)=X_ST(:,3)/10000;
scatter(X_ST(:,5)/5000,Res/(2*10^6),100,C,"filled"); 
set(ax1,'FontName','Palatino Linotype','FontSize',19,'FontWeight','normal','LineWidth',1);
xlabel('Killer T cells','FontSize',20,'FontWeight','normal','FontName','Palatino Linotype')
ylabel('Resource','FontSize',20,'FontWeight','normal','FontName','Palatino Linotype')
ylim([0 1])
C1=zeros(length(Res),3);
C1(:,1)=0;
C1(:,2)=(X_ST_post(:,1)+X_ST_post(:,2))/20000;
C1(:,3)=X_ST_post(:,3)/10000;
ylim([0 1])
figure;
ax2=axes();
scatter(X_ST_post(:,6)/5000,Res/(2*10^6),100,C1,"filled"); 
xlabel('Killer T cells','FontSize',20,'FontWeight','normal','FontName','Palatino Linotype')
ylabel('Resource','FontSize',20,'FontWeight','normal','FontName','Palatino Linotype')
set(ax2,'FontName','Palatino Linotype','FontSize',19,'FontWeight','normal','LineWidth',1);
ylim([0 1])


C=zeros(length(Res),3);
C=Res/max(Res);                                                                                        % The redness of each point denotes the associated resource with holding capacity
subplot(1,2,1)
scatter3((X_ST(:,1)+X_ST(:,2))/20000,X_ST(:,3)/10000,X_ST(:,5)/5000,100,C,"filled");                        % Plotting the tumor cell states and killer T cell components of the pre-ICI attractors
xlim([0 max((X_ST(:,1)+X_ST(:,2))/20000)])
ylim([0 max((X_ST(:,3))/10000)])
xlabel('PDL1- tumor cells')
ylabel('PDL1+ tumor cells')
zlabel('Killer T cells')
subplot(1,2,2)
scatter3((X_ST_post(:,1)+X_ST_post(:,2))/20000,X_ST_post(:,3)/10000,X_ST_post(:,6)/5000,100,C,"filled");   % Plotting the tumor cell states and killer T cell components of the pre-ICI attractors
xlim([0 max((X_ST_post(:,1)+X_ST_post(:,2))/20000)])
ylim([0 max((X_ST_post(:,3))/10000)])
zlim([0 max(X_ST_post(:,6)/5000)])
xlabel('PDL1- tumor cells')
ylabel('PDL1+ tumor cells')
zlabel('Killer T cells')