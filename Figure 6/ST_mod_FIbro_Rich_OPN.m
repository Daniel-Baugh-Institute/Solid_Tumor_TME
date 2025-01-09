function dydt=ST_mod_FIbro_Rich_OPN(t,y,u,P,K_knock)

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

K_CAFCR=P(11);             % CAF-driven growth of the tumor non-T exposed cells.

K_TKC=P(12);               % T cells-driven apoptosis of tumor cells.

K_TKCIFNG=P(13);           % IFNG induced T cells-driven apoptosis of tumor cells.

alpha=P(14);               % CAF barrier

K_CAFB=P(15);              % Barrier formation rate

K_CPDNPD= P(16);           % NPDL1 to PDL1 conversion

K_CNPDL1D=P(17);           % Death rate of NPDL1

K_CPDL1D=P(18);            % Death rate of PDL1

K_ResD=P(19);              % Degradation of Resources

Delta=P(20);               % Width of the barrier

%% T cell proliferation, interactions, and death

Y_TKNPD1=P(21);            % Carrying capacity for PD1+ T cells

K_TKPD=P(22);              % Proliferation rate of PD1+ T cells 

K_TKNPD=P(23);             % Proliferation rate of PD1- T cells

K_TKPDNPD1=P(24);          % Conversion from PD1+ to PD1- T cells

K_TKPDTEX=P(25);           % Conversion from PD1+ T cels to exhausted T cells  

K_CNPDT=P(26);             % MHC sensing of tumor cells for immune activation

K_CAFT=P(27);              % CAF-driven inhibition of T cells through regulatory T cells.

K_TKPDD=P(28);             % Death rate for killer PD1+ T cells

K_TEXD=P(29);              % Death rate for exhausted T cells

%% CAF proliferation, interactions, and death

Y_CAF=P(30);               % maximum carrying capacity of Fibroblasts

K_CAF=P(31);               % Proliferation rate of CAF cells 

K_OPNCAF=P(32);            % OPN-induced Proliferation rate of CAF cells

K_CTCAF=P(33);             % Tumor cells-induced growth of CAF

K_CTCAFR=P(34);            % Proximity factor for CAF and resistant tumor cells  

K_M2CAF=P(35);             % M2-driven proliferation rate for CAFs 

K_CAFD=P(36);              % Death rate for CAF

%% Macrophages
Y_M= P(37);               % Carrying capacity

K_M1= P(38);               % Growth rate of M1 phase macrophage

K_CANM1= P(39);            % Tumor cells-driven proliferation of macrophages

K_M1M2=P(40);              % Change of polarity

K_M1D= P(41);              % Death rate of M1 phase macrophage 

K_M2= P(42);               % Growth rate of M2 phase macrophage

K_CAFM2= P(43);            % CAF-driven proliferation of M2 macrophages

alpha_TM= P(44);           % Functional approximation of ICAM1 wrt T cells

K_M2D= P(45);              % Death rate of M2 phase macrophage 

%% IL-8 Secretion
K_M2IL8=P(46);             % IL8 secretion by M2

K_CAFIL8=P(47);            % IL8 secretion by CAF

K_CANIL8=P(48);            % IL8 secretion by tumor cells

K_IL8D=P(49);              % Degradation rate of IL8

%% OPN Secretion

K_CAFOPN=P(50);            % OPN secretion by CAF

K_CANOPN=P(51);            % OPN secretion by tumor cells

K_OPND=P(52);              % Degradation rate of OPN

K_OPNknock=P(53);          % Knock_out reaction rate of OPN
%% IFNG Secretion

K_TIFNG=P(54);             % IFNG secretion by T

K_IFNGD=P(55);             % Degradation rate of IFNG


%% Cell states, cyto and chemokines

CNPDL1=y(1);               % T_exposed Tumor stem cells

CPDL1=y(2);                % Tumor cells exposed to T-cells with pdl1

CRNPDL1=y(3);              % Tumor cells Non exposed to T-cells without pdl1

CRPDL1=y(4);               % Tumor cells Non exposed to T-cells with pdl1

Res=y(5);                  % Resource concentration

TKPD1= y(6);               % Killer T cells with pd1

TKNPD1= y(7);              % Killer T cells without pd1

TEX= y(8);                 % Exhausted T cells

CAF= y(9);                 % Cancer associated fibroblasts

M1=y(10);                  % M1 phase-macrophage

M2=y(11);                  % M2 Phase macrophage

IL8=y(12);                 % Interleukin 8

OPN=y(13);                 % Osteopontin

IFNG=y(14);                % Interferon gamma
%% Hyper parameters
I=tanh(alpha*K_CAFB*CAF);

alpha_T=exp(-Delta*alpha^2*CAF^2);

%% Population of cancer cell states

% PDL1-ve T-exposed tumor cell

F_ResCNPDL1=Res*K_RCNPDL1*CNPDL1*(1-CNPDL1/(Y_CNPDL1M*(1-I)+1))*1/(CPDL1+CRPDL1+CRNPDL1+1);                                  % Resource-driven proliferation of tumor stem cells

F_TxCNPDL1=(1+K_TXC*TEX/(1+TEX));                                                                                          % Exhausted T-cells-driven proliferation modulator of tumor stem cells 

F_CAFCNPDL1=(1+K_CAFC*CAF/(CAF+1));                                                                                        % CAF-driven proliferation of the tumor cells

F_CNPDL1CPDL1= K_CPDNPD*(IFNG/(IFNG+1))*CNPDL1;                                                                                            % Conversion from stem to PDL1, T-exposed tumor cells
 
F_TKCNPDL1= K_TKC*(1+K_TKCIFNG*IFNG/(IFNG+1))*CNPDL1*(TKPD1+TKNPD1);                                                                                   % T-cell driven apoptosis of the stem cell

F_DCNPDL1= (K_CNPDL1D+M1/(M1+1)+IFNG/(IFNG+1))*CNPDL1;                                                                                               % Death of stem cells

dydt_CNPDL1=F_ResCNPDL1*F_TxCNPDL1*F_CAFCNPDL1-F_TKCNPDL1-F_CNPDL1CPDL1-F_DCNPDL1;


% PDL1+, T-exposed tumor cells
    
F_ResCPDL1=Res*K_RCPDL1*CPDL1*(1-CPDL1/(Y_CPDL1M*(1-I)+1))*1/(CNPDL1+CRNPDL1+CRPDL1+1);                                     % Resource-driven proliferation of PDL1 +ve cells

F_TxCNPDL1=(1+K_TXC*TEX/(1+TEX));                                                                                         % Exhausted T-cells-driven proliferation modulator of tumor stem cells

F_CAFCNPDL1=(1+K_CAFC*CAF/(CAF+1));                                                                                       % CAF-driven proliferation of the tumor cells
   
F_TKNPDCPDL1=K_TKC*(1+K_TKCIFNG*IFNG/(IFNG+1))*CPDL1*TKNPD1;                                                                                          % Killer, PD1-ve T-cell driven death of PDL1+ve tumor cells

F_DCPDL1= (K_CPDL1D+M1/(M1+1)+IFNG/(IFNG+1))*CPDL1;                                                                                                 % Death of PDL1, T-exposed tumor cells   

dydt_CPDL1=F_ResCPDL1*F_TxCNPDL1*F_CAFCNPDL1+F_CNPDL1CPDL1-F_TKNPDCPDL1-F_DCPDL1;

% PDL1-ve Non-T-exposed tumor stem cell
  
F_ResCRNPDL1=Res*K_RCRNPDL1*CRNPDL1*(1-CRNPDL1/(Y_CNPDL1M*I+1))*1/(CPDL1+CRPDL1+CNPDL1+1);                                 % Resource-driven proliferation of tumor stem cells

F_CAFCRNPDL1=(1+K_CAFCR*CAF/(CAF+1));                                                                                    % CAF-driven proliferation of the tumor cells

F_CRNPDL1CRPDL1= K_CPDNPD*(IFNG/(IFNG+1))*CRNPDL1;                                                                                       % Conversion from stem to PDL1, T-exposed tumor cells

F_TxCR=(1+K_TXC*TEX/(1+TEX)*alpha_T);                                                                                    % Exhausted T-cells-driven proliferation modulator of tumor stem cells

F_TKCRNPDL1= K_TKC*(1+K_TKCIFNG*IFNG/(IFNG+1))*CRNPDL1*(TKPD1+TKNPD1)*alpha_T;                                                                       % T-cell driven apoptosis of the stem cell

F_DCRNPDL1= (K_CNPDL1D+M1/(M1+1)+IFNG/(IFNG+1))*CRNPDL1;                                                                                           % Death of stem cells

dydt_CRNPDL1=F_ResCRNPDL1*F_CAFCRNPDL1*F_TxCR-F_TKCRNPDL1-F_CRNPDL1CRPDL1-F_DCRNPDL1;

% PDL1+ve Non-T-exposed tumor stem cell

F_ResCRPDL1=Res*K_RCRPDL1*CRPDL1*(1-CRPDL1/(Y_CPDL1M*I+1))*1/(CPDL1+CRNPDL1+CNPDL1+1);                                    % Resource-driven proliferation of tumor stem cells

F_TKCRPDL1= K_TKC*(1+K_TKCIFNG*IFNG/(IFNG+1))*CRPDL1*(TKNPD1)*alpha_T;                                                                              % T-cell driven apoptosis of the stem cell

F_DCRPDL1= (K_CNPDL1D+M1/(M1+1)+IFNG/(IFNG+1))*CRPDL1;                                                                                            % Death of stem cells

dydt_CRPDL1=F_ResCRPDL1*F_CAFCRNPDL1*F_TxCR-F_TKCRPDL1+F_CRNPDL1CRPDL1-F_DCRPDL1;

%% Resource concentration

F_ProRes=K_RIn*(Y_RM-Res);

F_ConsRes=Res*(K_RCNPDL1*(CNPDL1+CRNPDL1)+K_RCPDL1*(CPDL1+CRPDL1));

F_DRes= K_ResD*Res;

dydt_Res=F_ProRes-F_ConsRes-F_DRes;

%% T-cell population

% Killer PD1+ T cells

F_ProTKPD1=K_TKPD*TKPD1*(1-TKPD1/Y_TKNPD1);                                                                                                 % Proliferation of killer, PD1+ve T-cells

F_TKPD1TKNPD1=K_TKPDNPD1*u*TKPD1;                                                                                                           % ICI-based conversion from PD1+ve to -ve killer T cells

F_TKPD1CPDL1=(1+K_CNPDT*(CNPDL1+CRNPDL1*alpha_T)/(CNPDL1+CRNPDL1*alpha_T+1)*K_CAFT/(1+CAF));                                                % Tumor cells-driven activation of T cells

F_TKPD1TEX=K_TKPDTEX*TKPD1*(M2*(CPDL1+CRPDL1*alpha_T))/(M2*(CPDL1+CRPDL1*alpha_T)+1);                                                       % PDL1+ve and M2 driven conversion of PD1+ve to exhausted T cells  

%F_TKCAF=(CAF+1)/((K_CAFT+1)*CAF+1);

F_DTKPD1= K_TKPDD*TKPD1;                                                                                                                    % Death rate of PD1+, killer T cells

dydt_TKPD1 = F_ProTKPD1*F_TKPD1CPDL1-F_TKPD1TKNPD1-F_TKPD1TEX-F_DTKPD1; 

% Killer PD1- T cells

F_ProTKNPD1=K_TKNPD*TKNPD1*(1-TKNPD1/Y_TKNPD1);                                                                                              % Proliferation of killer, PD1+ve T-cells
 
F_TKPD1CPDL1=(1+K_CNPDT*(CNPDL1+CPDL1+(CRNPDL1+CRPDL1)*alpha_T)/(CNPDL1+CPDL1+(CRNPDL1+CRPDL1)*alpha_T+1)*K_CAFT/(1+CAF));                  % Tumor cells-driven activation of T cells

F_DTKNPD1= K_TKPDD*TKNPD1;

%dydt_TKNPD1=F_TKPD1TKNPD1-F_DTKNPD1;

dydt_TKNPD1=F_ProTKNPD1*u*F_TKPD1CPDL1+F_TKPD1TKNPD1-F_DTKNPD1;

% Exhausted T cells
F_DTEX=K_TEXD*TEX;                                                                                                                          % Death rate of exhausted t cells

dydt_TEX=F_TKPD1TEX-F_DTEX;

%% CAF population

F_ProCAF= K_CAF*CAF*(1-CAF/Y_CAF);                                                                                                            % Proliferation of CAF

F_OPNCAF=(1+K_OPNCAF*OPN/(OPN+1));

F_CANCAF=(1+ K_CTCAF*(CPDL1+CNPDL1)/(10^4+CPDL1+CNPDL1));                                                                                     % ICI-based conversion from PD1+ve to -ve killer T cells

F_CANRCAF=(1+ K_CTCAFR*(CRPDL1+CRNPDL1)/(10^4+CRPDL1+CRNPDL1));                                                                                % ICI-based conversion from PD1+ve to -ve killer T cells

F_M2CAF=(1+K_M2CAF*M2/(M2+10^4));

F_DCAF = K_CAFD*CAF;                                                                                                                            % Death rate of CAFs

dydt_CAF = F_ProCAF*F_CANCAF*F_CANRCAF*F_M2CAF*F_OPNCAF-F_DCAF; 

%% Macrophages population
% M1 Phase
F_ProM1 = K_M1*M1*(1-M1/(Y_M-M2));                                                                                                         % Proliferation of CAF

F_CANM1 =(1+ K_CANM1*(CPDL1+CNPDL1+CRPDL1+CRNPDL1)/(1+CPDL1+CNPDL1+CRPDL1+CRNPDL1));                                                       % ICI-based conversion from PD1+ve to -ve killer T cells

F_M1M2=K_M1M2*(alpha_TM*(TKPD1+TKNPD1))^2/((alpha_TM*(TKPD1+TKNPD1))^2+1)*M1;

F_DM1 = K_M1D*M1;                                                                                                                          % Death rate of CAFs

dydt_M1 = F_ProM1*F_CANM1-F_M1M2-F_DM1;

% M2 Phase
F_ProM2 = K_M2*M2*(1-M2/(Y_M-M1));                                                                                                         % Proliferation of CAF

F_CAFM2 =(1+K_CAFM2*CAF/(CAF+1));

F_DM2 = K_M2D*M2;                                                                                                                          % Death rate of CAFs

dydt_M2 = F_ProM2*F_CAFM2+F_M1M2-F_DM2;

%% IL-8 concentration
F_M2IL8= K_M2IL8*M2;                                                                                                                         % IL8 secretion from M2

F_CAFIL8= K_CAFIL8*CAF;                                                                                                                     % IL8 secretion from CAF

F_CANIL8= K_CANIL8*(CPDL1+CNPDL1+CRPDL1+CRNPDL1);                                                                                            % I8 secretion from Tumor cells

F_DIL8= K_IL8D*IL8;

dydt_IL8= F_M2IL8+F_CAFIL8+F_CANIL8-F_DIL8;

%% OPN concentration
F_CAFOPN= K_CAFOPN*CAF;                                                                                                                         % OPN secretion from CAF

F_CANOPN= K_CANOPN*(CPDL1+CNPDL1+CRPDL1+CRNPDL1);                                                                                               % OPN secretion from Tumor cells

F_DOPN= K_OPND*OPN;

F_PhiOPN= K_knock*OPN*K_OPNknock;

dydt_OPN= F_CAFOPN+F_CANOPN-F_DOPN-F_PhiOPN;

%% IFNG concentration
F_TIFNG= K_TIFNG*(TKNPD1+TKPD1);                                                                                                                         % OPN secretion from CAF

F_OPNIFNG= 1/(OPN+1);                                                                                                                         % OPN secretion from CAF

F_DIFNG= K_IFNGD*IFNG;

dydt_IFNG= F_TIFNG*F_OPNIFNG-F_DIFNG;

dydt = [dydt_CNPDL1;dydt_CPDL1;dydt_CRNPDL1;dydt_CRPDL1;dydt_Res;dydt_TKPD1;dydt_TKNPD1;dydt_TEX;dydt_CAF;dydt_M1;dydt_M2;dydt_IL8;dydt_OPN;dydt_IFNG];

end