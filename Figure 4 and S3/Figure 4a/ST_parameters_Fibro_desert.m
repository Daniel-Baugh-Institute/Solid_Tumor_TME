function P=ST_parameters_Fibro_desert

%% Resource to cancer cells
K_RCST=100;         %  per unit of resource per unit of time. Resource to tumor stem cell
P(1)=K_RCST;       

K_RCNPD=80;         %  per unit of resource per unit of time. Resource to PDL1- tumor cell
P(2)=K_RCNPD;

K_RCPD=10;          %  per unit of resource per unit of time. Resource to PDL1+ tumor cell
P(3)=K_RCPD;

K_RIn=1000;         %  Resource per unit time. Default resource intake rate  
P(4)=K_RIn;

%% Cancer cell proliferation, interactions, and death
Y_CSTM=10^4;        % Cells. Carrying capacity of stem tumor cells  
P(5)=Y_CSTM;

Y_CNPDM=10^4;       % Cells. Carrying capacity of PDL1- tumor cells
P(6)=Y_CNPDM;

Y_CPDM=10^4;        % Cells. Carrying capacity of PDL1+ tumor cells
P(7)=Y_CPDM;

Y_RM=10^4;          % Resources. Maximum resource withholding capacity.
P(8)=Y_RM;

K_TXC=20;           % Exhausted T cells-driven growth of tumor cells.
P(9)=K_TXC;

K_TKC=0.3;          % T cells-driven apoptosis of tumor cells.
P(10)=K_TKC;

K_CSTNPD=5;         % Stem to NPDL1 conversion.
P(11)=K_CSTNPD;

K_CSTPD=7;          % Stem to PDL1 conversion.
P(12)=K_CSTPD;

K_CPDNPD=2;         % NPDL1 to PDL1 conversion
P(13)=K_CPDNPD;

K_CSTD=0.002;           % Death rate of Stem
P(14)=K_CSTD;

alpha_C=1;           % Spatial competition 
P(15)=alpha_C;

K_CNPD=6;           % Death rate of NPDL1
P(16)=K_CNPD;

K_CPD=6;            % Death rate of PDL1
P(17)=K_CPD;

K_ResD=8;           % Degradation of Resources
P(18)=K_ResD;

Y_TKPD1=5000;       % Carrying capacity for PD1+ T cells
P(19)=Y_TKPD1;

Y_TKNPD1=5000;      % Carrying capacity for PD1- T cells
P(20)=Y_TKNPD1;

%% T cell proliferation, interactions, and death
K_TKPD=40;          % Proliferation rate of PD1+ T cells 
P(21)=K_TKPD;

K_TKNPD=45;         % Proliferation rate of PD1- T cells
P(22)=K_TKNPD;

K_TKPDNPD1=1.5*10^4;     % Conversion from PD1+ to PD1- T cells
P(23)=K_TKPDNPD1;

K_TKPDTEX=8;       % Conversion from PD1+ T cels to exhausted T cells  
P(24)=K_TKPDTEX;

K_CNPDT=40;         % MHC sensing of tumor cells for immune activation
P(25)=K_CNPDT;

K_TKPDD=8;          % Death rate for killer PD1+ T cells
P(26)=K_TKPDD;

K_TKNPDD=8;         % Death rate for killer PD1+ T cells
P(27)=K_TKNPDD;

K_TEXD=10;           % Death rate for exhausted T cells
P(28)=K_TEXD;

end