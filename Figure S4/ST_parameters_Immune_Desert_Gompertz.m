function P=ST_parameters_Immune_Desert_Gompertz

%% Load the parameters

%% Resource to cancer cells

K_RCNPDL1=100;
P(1)=K_RCNPDL1;            %  per unit of resource per unit of time. Resource to T-exposed PDL1- tumor cell 

K_RCRNPDL1=100;
P(2)=K_RCRNPDL1;           %  per unit of resource per unit of time. Resource to Non-T-exposed PDL1- tumor cell

K_RCPDL1=100;
P(3)=K_RCPDL1;             %  per unit of resource per unit of time. Resource to T-exposed PDL1+ tumor cell

K_RCRPDL1=10;
P(4)=K_RCRPDL1;             %  per unit of resource per unit of time. Resource to Non-T-exposed PDL1+ tumor cell

K_RIn=2000;
P(5)=K_RIn;                 %  Resource per unit time. Default resource intake rate  

%% Cancer cell proliferation, interactions, and death
Y_CNPDL1M=10^4;
P(6)=Y_CNPDL1M;             % Carrying capacity for PDL1- cells

Y_CPDL1M=10^4;
P(7)=Y_CPDL1M;             % Carrying capacity for PDL1+ T cells

Y_RM=10^4;
P(8)=Y_RM;                 % Resources. Maximum resource withholding capacity.

K_TXC=20;
P(9)=K_TXC;                % Exhausted T cells-driven growth of tumor cells.

K_CAFC=60;
P(10)=K_CAFC;             % CAF-driven growth of the tumor cells.

K_TKC=50;
P(11)=K_TKC;               % T cells-driven apoptosis of tumor cells.

K_CPDNPD=0.25;
P(12)=K_CPDNPD;           % NPDL1 to PDL1 conversion

K_CNPDL1D=12;
P(13)=K_CNPDL1D;           % Death rate of NPDL1

K_CPDL1D=12;
P(14)=K_CPDL1D;            % Death rate of PDL1

K_ResD=20;
P(15)=K_ResD;              % Degradation of Resources

%% T cell proliferation, interactions, and death

Y_TKNPD1=5000;
P(16)=Y_TKNPD1;            % Carrying capacity for PD1+ T cells

K_TKPD=5;
P(17)=K_TKPD;              % Proliferation rate of PD1+ T cells 

K_TKNPD=10;
P(18)=K_TKNPD;             % Proliferation rate of PD1- T cells

K_TEX=5;                   
P(19)=K_TEX;               % Proliferation of exhausted T cells 

K_TKPDNPD1=500;
P(20)=K_TKPDNPD1;          % Conversion from PD1+ to PD1- T cells

K_TKPDTEX=20;
P(21)=K_TKPDTEX;           % Conversion from PD1+ T cels to exhausted T cells  

K_CNPDT=0.2;
P(22)=K_CNPDT;             % MHC sensing of tumor cells for immune activation

K_IL2TK=0.25;                
P(23)=K_IL2TK;             % IL2-induced proliferation of killer T cells

K_CAFT=0.2;
P(24)=K_CAFT;              % CAF-driven inhibition of T cells through regulatory T cells.

K_TKPDD=36;
P(25)=K_TKPDD;             % Death rate for killer PD1+ T cells

K_TEXD=25;
P(26)=K_TEXD;              % Death rate for exhausted T cells

%% CAF proliferation, interactions, and death
Y_CAF=5000;
P(27)=Y_CAF;               % maximum carrying capacity of Fibroblasts

K_CAF=50;
P(28)=K_CAF;               % Proliferation rate of CAF cells 

K_CTCAF=50;
P(29)=K_CTCAF;             % Tumor cells-induced growth of CAF

K_M2CAF=10;
P(30)=K_M2CAF;             % M2-driven proliferation rate for CAFs

K_CAFD=7;
P(31)=K_CAFD;              % Death rate for CAF

%% Macrophage proliferation, interactions, and death
Y_M=5000;
P(32)=Y_M;               % Carrying capacity

K_M1=10;
P(33)=K_M1;               % Growth rate of M1 phase macrophage

K_CANM1=10;
P(34)=K_CANM1;            % Tumor cells-driven proliferation of macrophages

K_M1M2=40;
P(35)=K_M1M2;              % Change of polarity

K_M1D=8;
P(36)=K_M1D;              % Death rate of M1 phase macrophage 

K_M2=10;
P(37)=K_M2;               % Growth rate of M2 phase macrophage

K_CAFM2=150;
P(38)=K_CAFM2;            % CAF-driven proliferation of M2 macrophages

alpha_TM=0.0005;          
P(39)=alpha_TM;           % Functional approximation of ICAM1 wrt T cells

K_M2D=8;
P(40)=K_M2D;              % Death rate of M2 phase macrophage 

%% IL-2 Concentration
K_TKIL2=0.5;
P(41)=K_TKIL2;             % IL8 secretion by M2

K_IL2D=2;                  % Degradation rate of IL8
P(42)=K_IL2D;
end