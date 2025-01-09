function dydt=ST_mod_Immune_Rich(t,y,u,u_Tk,P)

%% Load the parameters

%% Resource to cancer cells

K_RCNPDL1=P(1);            %  per unit of resource per unit of time. Resource to T-exposed PDL1- tumor cell 

K_RCRNPDL1=P(2);           %  per unit of resource per unit of time. Resource to Non-T-exposed PDL1- tumor cell

K_RCPDL1=P(3);             %  per unit of resource per unit of time. Resource to T-exposed PDL1+ tumor cell

K_RCRPDL1=P(4);            %  per unit of resource per unit of time. Resource to Non-T-exposed PDL1+ tumor cell

K_RIn=P(5);                %  Resource per unit time. Default resource intake rate  

%% Cancer cell proliferation, interactions, and death

Y_CNPDL1M=P(6);            % Carrying capacity for PD1- T cells

Y_CPDL1M=P(7);             % Carrying capacity for PD1+ T cells

Y_RM=P(8);                 % Resources. Maximum resource withholding capacity.

K_TXC=P(9);                % Exhausted T cells-driven growth of tumor cells.

K_CAFC= P(10);             % CAF-driven growth of the tumor cells.

K_TKC=P(11);               % T cells-driven apoptosis of tumor cells.

K_CPDNPD= P(12);           % NPDL1 to PDL1 conversion

K_CNPDL1D=P(13);           % Death rate of NPDL1

K_CPDL1D=P(14);            % Death rate of PDL1

K_ResD=P(15);              % Degradation of Resources

%% T cell proliferation, interactions, and death

Y_TKNPD1=P(16);            % Carrying capacity for PD1+ T cells

K_TKPD=P(17);              % Proliferation rate of PD1+ T cells 

K_TKNPD=P(18);             % Proliferation rate of PD1- T cells

K_TKPDNPD1=P(19);          % Conversion from PD1+ to PD1- T cells

K_TKPDTEX=P(20);           % Conversion from PD1+ T cels to exhausted T cells  

K_CNPDT=P(21);             % MHC sensing of tumor cells for immune activation

K_IL2TK=P(22);             % IL2-induced proliferation of killer T cells

K_CAFT=P(23);              % CAF-driven inhibition of T cells through regulatory T cells.

K_TKPDD=P(24);             % Death rate for killer PD1+ T cells

K_TEXD=P(25);              % Death rate for exhausted T cells

%% CAF proliferation, interactions, and death

Y_CAF=P(26);               % maximum carrying capacity of Fibroblasts

K_FWT=P(27);               % Proliferation rate of wild type fibroblasts

K_FWTCAF=P(28);            % FWT to CAF conversion

K_CAFFWT=P(29);            % CAF to FWT conversion

alpha_FWTC=P(30);          % Proximity of FWT and Tumor cells

K_FWTD=P(31);              % Death rate of FWT

K_CAF=P(32);               % Proliferation rate of CAF cells 

K_CTCAF=P(33);             % Tumor cells-induced growth of CAF

K_M2CAF=P(34);             % M2-driven proliferation rate for CAFs 

K_CAFD=P(35);              % Death rate for CAF

%% Macrophages
Y_M= P(36);               % Carrying capacity

K_M1= P(37);               % Growth rate of M1 phase macrophage

K_CANM1= P(38);            % Tumor cells-driven proliferation of macrophages

K_M1M2=P(39);              % Change of polarity

K_M1D= P(40);              % Death rate of M1 phase macrophage 

K_M2= P(41);               % Growth rate of M2 phase macrophage

K_CAFM2= P(42);            % CAF-driven proliferation of M2 macrophages

alpha_TM= P(43);           % Functional approximation of ICAM1 wrt T cells

K_M2D= P(44);              % Death rate of M2 phase macrophage 

%% IL8 secretion

K_M2IL8=P(45);             % IL8 secretion by M2

K_CAFIL8=P(46);            % IL8 secretion by CAF

K_CANIL8=P(47);            % IL8 secretion by tumor cells

K_IL8D=P(48);              % Degradation rate of IL8

%% IL2 secretion

K_TKIL2=P(49);             % IL8 secretion by M2

K_IL2D=P(50);              % Degradation rate of IL8

%% Cell states, cyto and chemokines

CNPDL1=y(1);               % T_exposed Tumor stem cells

CPDL1=y(2);                % Tumor cells exposed to T-cells with pdl1

Res=y(3);                  % Resource concentration

TKPD1= y(4);               % Killer T cells with pd1

TKNPD1= y(5);              % Killer T cells without pd1

TEX= y(6);                 % Exhausted T cells

FWT=y(7);                  % Wild type fibroblasts

CAF= y(8);                 % Cancer associated fibroblasts

M1=y(9);                   % M1 phase-macrophage

M2=y(10);                   % M2 Phase macrophage

IL8=y(11);                 % IL8 Concentration

IL2=y(12);                 % IL-2 Concentration
%% Population of cancer cell states

% PDL1-ve T-exposed tumor cell

F_ResCNPDL1=Res*K_RCNPDL1*CNPDL1*(1-CNPDL1/(Y_CNPDL1M))*1/(CPDL1+1);                                  % Resource-driven proliferation of tumor stem cells

F_TxCNPDL1=(1+K_TXC*TEX/(1+TEX));                                                                                          % Exhausted T-cells-driven proliferation modulator of tumor stem cells 

F_CAFCNPDL1=(1+K_CAFC*CAF/(CAF+1));                                                                                        % CAF-driven proliferation of the tumor cells

F_CNPDL1CPDL1= K_CPDNPD*CNPDL1;                                                                                            % Conversion from stem to PDL1, T-exposed tumor cells
 
F_TKCNPDL1= K_TKC*CNPDL1*(TKPD1+TKNPD1);                                                                                   % T-cell driven apoptosis of the stem cell

F_DCNPDL1= (K_CNPDL1D+M1/(M1+1))*CNPDL1;                                                                                   % Death of stem cells

dydt_CNPDL1=F_ResCNPDL1*F_TxCNPDL1*F_CAFCNPDL1-F_TKCNPDL1-F_CNPDL1CPDL1-F_DCNPDL1;


% PDL1+, T-exposed tumor cells
    
F_ResCPDL1=Res*K_RCPDL1*CPDL1*(1-CPDL1/(Y_CPDL1M))*1/(CNPDL1+1);                                     % Resource-driven proliferation of PDL1 +ve cells

F_TxCNPDL1=(1+K_TXC*TEX/(1+TEX));                                                                                         % Exhausted T-cells-driven proliferation modulator of tumor stem cells

F_CAFCNPDL1=(1+K_CAFC*CAF/(CAF+1));                                                                                       % CAF-driven proliferation of the tumor cells
   
F_TKNPDCPDL1=K_TKC*CPDL1*TKNPD1;                                                                                          % Killer, PD1-ve T-cell driven death of PDL1+ve tumor cells

F_DCPDL1= (K_CPDL1D+M1/(M1+1))*CPDL1;                                                                                                 % Death of PDL1, T-exposed tumor cells   

dydt_CPDL1=F_ResCPDL1*F_TxCNPDL1*F_CAFCNPDL1+F_CNPDL1CPDL1-F_TKNPDCPDL1-F_DCPDL1;

%% Resource concentration

F_ProRes=K_RIn*(Y_RM-Res);

F_ConsRes=Res*(K_RCNPDL1*(CNPDL1+CPDL1));
F_DRes= K_ResD*Res;

dydt_Res=F_ProRes-F_ConsRes-F_DRes;

%% T-cell population

% Killer PD1+ T cells

F_ProTKPD1=K_TKPD*TKPD1*(1-TKPD1/Y_TKNPD1);                                                                                                 % Proliferation of killer, PD1+ve T-cells

F_TKPD1TKNPD1=K_TKPDNPD1*u*TKPD1;                                                                                                           % ICI-based conversion from PD1+ve to -ve killer T cells

F_TKPD1CPDL1=(1+K_CNPDT*(CNPDL1)/(CNPDL1+1)*K_CAFT/(1+CAF));                                                % Tumor cells-driven activation of T cells

F_TKPD1TEX=K_TKPDTEX*TKPD1*(M2*(CPDL1))/(M2*(CPDL1)+1);                                                       % PDL1+ve and M2 driven conversion of PD1+ve to exhausted T cells  

F_TKIL2=(1+K_IL2TK*(IL2)/(IL2+1));
%F_TKCAF=(CAF+1)/((K_CAFT+1)*CAF+1);

F_DTKPD1= K_TKPDD*TKPD1;                                                                                                                    % Death rate of PD1+, killer T cells

dydt_TKPD1 = F_ProTKPD1*F_TKPD1CPDL1*F_TKIL2-F_TKPD1TKNPD1-F_TKPD1TEX-F_DTKPD1+u_Tk; 

% Killer PD1- T cells

F_ProTKNPD1=K_TKNPD*TKNPD1*(1-TKNPD1/Y_TKNPD1);                                                                                              % Proliferation of killer, PD1+ve T-cells
 
F_TKPD1CPDL1=(1+K_CNPDT*(CNPDL1+CPDL1)/(CNPDL1+CPDL1+1)*K_CAFT/(1+CAF));                                                                     % Tumor cells-driven activation of T cells

F_DTKNPD1= K_TKPDD*TKNPD1;

%dydt_TKNPD1=F_TKPD1TKNPD1-F_DTKNPD1;

dydt_TKNPD1=F_ProTKNPD1*u*F_TKPD1CPDL1*F_TKIL2+F_TKPD1TKNPD1-F_DTKNPD1;

% Exhausted T cells
F_DTEX=K_TEXD*TEX;                                                                                                                          % Death rate of exhausted t cells

dydt_TEX=F_TKPD1TEX-F_DTEX;

%% Fibroblasts population
% Wild type fibroblasts

F_ProFWT=K_FWT*FWT*(1-FWT/(Y_CAF-CAF));                              % Proliferation of wild type fibroblasts

F_CAFFWT= K_CAFFWT*CAF;

F_FWTCAF= K_FWTCAF*(alpha_FWTC*(CPDL1+CNPDL1)/(CPDL1+CNPDL1+1))*FWT;

F_DFWT= K_FWTD*FWT;

dydt_FWT= F_ProFWT+F_CAFFWT-F_FWTCAF-F_DFWT;

% Cancer-associated fibroblasts

F_ProCAF= K_CAF*CAF*(1-CAF/(Y_CAF-FWT));                                                                                                         % Proliferation of CAF

F_CANCAF=(1+ K_CTCAF*(CPDL1+CNPDL1)/(10^4+(CPDL1+CNPDL1)));                                                                                     % ICI-based conversion from PD1+ve to -ve killer T cells

F_M2CAF=(1+K_M2CAF*M2/(M2+1));

F_DCAF = K_CAFD*CAF;                                                                                                                        % Death rate of CAFs

dydt_CAF = F_ProCAF*F_CANCAF*F_M2CAF+F_FWTCAF-F_CAFFWT-F_DCAF; 

%% Macrophages population
% M1 Phase
F_ProM1 = K_M1*M1*(1-M1/(Y_M-M2));                                                                                                         % Proliferation of CAF

F_CANM1 =(1+ K_CANM1*(CPDL1+CNPDL1)/(1+CPDL1+CNPDL1));                                                       % ICI-based conversion from PD1+ve to -ve killer T cells

F_M1M2=K_M1M2*(alpha_TM*(TKPD1+TKNPD1))^2/((alpha_TM*(TKPD1+TKNPD1))^2+1)*M1;

F_DM1 = K_M1D*M1;                                                                                                                          % Death rate of CAFs

dydt_M1 = F_ProM1*F_CANM1-F_M1M2-F_DM1;

% M2 Phase
F_ProM2 = K_M2*M2*(1-M2/(Y_M-M1));                                                                                                         % Proliferation of CAF

F_CAFM2 =(1+K_CAFM2*CAF/(CAF+10^5));

F_DM2 = K_M2D*M2;                                                                                                                          % Death rate of CAFs

dydt_M2 = F_ProM2*F_CAFM2+F_M1M2-F_DM2;

%% IL-8 concentration
F_M2IL8= K_M2IL8*M2;                                                                                                                         % IL8 secretion from M2

F_CAFIL8= K_CAFIL8*CAF;                                                                                                                     % IL8 secretion from CAF

F_CANIL8= K_CANIL8*(CPDL1+CNPDL1);                                                                                                          % I8 secretion from Tumor cells

F_DIL8= K_IL8D*IL8;

dydt_IL8= F_M2IL8+F_CAFIL8+F_CANIL8-F_DIL8;

%% IL-2 concentration
F_TKIL2= K_TKIL2*(TKNPD1+TKPD1);                                                                                                                         % IL8 secretion from M2

F_DIL2= K_IL2D*IL2;

dydt_IL2= F_TKIL2-F_DIL2;


dydt = [dydt_CNPDL1;dydt_CPDL1;dydt_Res;dydt_TKPD1;dydt_TKNPD1;dydt_TEX;dydt_FWT;dydt_CAF;dydt_M1;dydt_M2;dydt_IL8; dydt_IL2];

end