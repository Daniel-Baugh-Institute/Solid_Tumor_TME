function P=ST_parameters_Fibro_Rich_OPN(alpha)

%% Load the parameters

%% Resource to cancer cells

K_RCNPDL1=100;
P(1)=K_RCNPDL1;            %  per unit of resource per unit of time. Resource to T-exposed PDL1- tumor cell 

K_RCRNPDL1=100;
P(2)=K_RCRNPDL1;           %  per unit of resource per unit of time. Resource to Non-T-exposed PDL1- tumor cell

K_RCPDL1=10;
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

K_CAFC=50;
P(10)=K_CAFC;             % CAF-driven growth of the tumor cells.

K_CAFCR=50;
P(11)=K_CAFCR;             % CAF-driven growth of the tumor non-T exposed cells.

K_TKC=50;
P(12)=K_TKC;               % T cells-driven apoptosis of tumor cells.

K_TKCIFNG=3;             
P(13)=K_TKCIFNG;          % IFNG induced T cells-driven apoptosis of tumor cells.

P(14)=alpha;               % CAF barrier

alpha_C=1;                 % Competition index
P(15)=alpha_C;

K_CAFB=0.001;
P(16)=K_CAFB;              % Barrier formation rate

K_CPDNPD=8;
P(17)=K_CPDNPD;           % NPDL1 to PDL1 conversion

K_CNPDL1D=8;
P(18)=K_CNPDL1D;           % Death rate of NPDL1

K_CPDL1D=9;
P(19)=K_CPDL1D;            % Death rate of PDL1

K_ResD=10;
P(20)=K_ResD;              % Degradation of Resources

Delta=10^-3; 
P(21)=Delta;               % CAF barrier Thickness

%% T cell proliferation, interactions, and death
Y_TKNPD1=5000;
P(22)=Y_TKNPD1;            % Carrying capacity for PD1+ T cells

K_TKPD=45;
P(23)=K_TKPD;              % Proliferation rate of PD1+ T cells 

K_TKNPD=55;
P(24)=K_TKNPD;             % Proliferation rate of PD1- T cells

K_TKPDNPD1=500;
P(25)=K_TKPDNPD1;          % Conversion from PD1+ to PD1- T cells

K_TKPDTEX=12;
P(26)=K_TKPDTEX;           % Conversion from PD1+ T cels to exhausted T cells  

K_CNPDT=90;
P(27)=K_CNPDT;             % MHC sensing of tumor cells for immune activation

K_CAFT=0.2;
P(28)=K_CAFT;              % CAF-driven inhibition of T cells through regulatory T cells.

K_TKPDD=7;
P(29)=K_TKPDD;             % Death rate for killer PD1+ T cells

K_TEXD=9;
P(30)=K_TEXD;              % Death rate for exhausted T cells

%% CAF proliferation, interactions, and death
Y_CAF=5000;
P(31)=Y_CAF;               % maximum carrying capacity of Fibroblasts

K_CAF=1;
P(32)=K_CAF;               % Proliferation rate of CAF cells 

K_OPNCAF=50;
P(33)=K_OPNCAF;            % OPN-induced Proliferation rate of CAF cells

K_CTCAF=6;
P(34)=K_CTCAF;             % Tumor cells-induced growth of CAF

K_CTCAFR=0.02;
P(35)=K_CTCAFR;            % Proximity factor for CAF and resistant tumor cells  

K_M2CAF=1;
P(36)=K_M2CAF;             % M2-driven proliferation rate for CAFs

K_CAFD=12;
P(37)=K_CAFD;              % Death rate for CAF

%% Macrophage proliferation, interactions, and death
Y_M=5000;
P(38)=Y_M;                 % Carrying capacity

K_M1=10;
P(39)=K_M1;                % Growth rate of M1 phase macrophage

K_CANM1=10;
P(40)=K_CANM1;             % Tumor cells-driven proliferation of macrophages

K_M1M2=40;
P(41)=K_M1M2;              % Change of polarity

K_M1D=8;
P(42)=K_M1D;               % Death rate of M1 phase macrophage 

K_M2=10;
P(43)=K_M2;                % Growth rate of M2 phase macrophage

K_CAFM2=120;
P(44)=K_CAFM2;             % CAF-driven proliferation of M2 macrophages

alpha_TM=0.0005;          
P(45)=alpha_TM;            % Functional approximation of ICAM1 wrt T cells

K_M2D=8;
P(46)=K_M2D;               % Death rate of M2 phase macrophage 

%% IL-8 Concentration

K_M2IL8=10;            
P(47)=K_M2IL8;            % IL8 secretion by M2

K_CAFIL8=0.1;
P(48)=K_CAFIL8;           % IL8 secretion by CAF

K_CANIL8=0.1;
P(49)=K_CANIL8;           % IL8 secretion by tumor cells

K_IL8D=8;                  
P(50)=K_IL8D;             % Degradation rate of IL8
%% OPN Secretion

K_CAFOPN=2;
P(51)=K_CAFOPN;            % OPN secretion by CAF

K_CANOPN=0.05;
P(52)=K_CANOPN;            % OPN secretion by tumor cells

K_OPND=8;
P(53)=K_OPND;              % Degradation rate of OPN

K_OPNknock=30;
P(54)=K_OPNknock;          % Knock_out reaction rate of OPN
%% IFNG Secretion

K_TIFNG=5;
P(55)=K_TIFNG;            % IFNG secretion by T cells

K_IFNGD=8;
P(56)=K_IFNGD;            % Degradation rate of IFNG


end