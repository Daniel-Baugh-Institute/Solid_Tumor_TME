function dydt=ST_mod_FIbro_Desert(t,y,u,P)

%% Load the parameters

%% Resource to cancer cells

K_RCST=P(1);        %  per unit of resource per unit of time. Resource to tumor stem cell 

K_RCNPD=P(2);       %  per unit of resource per unit of time. Resource to PDL1- tumor cell

K_RCPD=P(3);        %  per unit of resource per unit of time. Resource to PDL1+ tumor cell

K_RIn=P(4);         %  Resource per unit time. Default resource intake rate  

%% Cancer cell proliferation, interactions, and death

Y_CSTM=P(5);        % Cells. Carrying capacity of stem tumor cells  

Y_CNPDM=P(6);       % Cells. Carrying capacity of PDL1- tumor cells

Y_CPDM=P(7);        % Cells. Carrying capacity of PDL1+ tumor cells

Y_RM=P(8);          % Resources. Maximum resource withholding capacity.

K_TXC=P(9);         % Exhausted T cells-driven growth of tumor cells.

K_TKC=P(10);        % T cells-driven apoptosis of tumor cells.

K_CSTNPD=P(11);     % Stem to NPDL1 conversion.

K_CSTPD=P(12);      % Stem to PDL1 conversion.

K_CPDNPD= P(13);    % NPDL1 to PDL1 conversion

K_CSTD=P(14);       % Death rate of Stem

alpha_C=P(15);      % Spatial competition

K_CNPD=P(16);       % Death rate of NPDL1

K_CPD=P(17);        % Death rate of PDL1

K_ResD=P(18);       % Degradation of Resources

Y_TKPD1=P(19);      % Carrying capacity for PD1+ T cells

Y_TKNPD1=P(20);     % Carrying capacity for PD1- T cells

%% T cell proliferation, interactions, and death
K_TKPD=P(21);       % Proliferation rate of PD1+ T cells 

K_TKNPD=P(22);      % Proliferation rate of PD1- T cells

K_TKPDNPD1=P(23);   % Conversion from PD1+ to PD1- T cells

K_TKPDTEX=P(24);    % Conversion from PD1+ T cels to exhausted T cells  

K_CNPDT=P(25);      % MHC sensing of tumor cells for immune activation

K_TKPDD=P(26);      % Death rate for killer PD1+ T cells

K_TKNPDD=P(27);     % Death rate for killer PD1+ T cells

K_TEXD=P(28);       % Death rate for exhausted T cells

%% cell states, cyto and chemokines
CST= y(1);          % T_exposed Tumor stem cells

CNPDL1=y(2);        % Tumor cells exposed to T-cells without pdl1   

CPDL1=y(3);         % Tumor cells exposed to T-cells with pdl1

Res=y(4);           % Resource concentration

TKPD1= y(5);        % Killer T cells with pd1

TKNPD1= y(6);       % Killer T cells without pd1

TEX= y(7);          % Exhausted T cells

%% Population of cancer cell states

% T-exposed tumor stem cell

F_ResCST=Res*K_RCST*CST*(1-CST/Y_CSTM)*1/(alpha_C*(CNPDL1+CPDL1)+1);                                   % Resource-driven proliferation of tumor stem cells

F_TxCST=(1+K_TXC*TEX/(1+TEX));                                                               % Exhausted T-cells-driven proliferation modulator of tumor stem cells 

F_CSTCNPDL1=K_CSTPD*CST;                                                                     % Conversion from stem to non-PDL1, T-exposed tumor cells

F_CSTCPDL1=K_CSTNPD*CST;                                                                     % Conversion from stem to PDL1, T-exposed tumor cells
 
F_TKCST= K_TKC*CST*(TKPD1+TKNPD1);                                                           % T-cell driven apoptosis of the stem cell

F_DCST= K_CSTD*CST;                                                                          % Death of stem cells

dydt_CST=F_ResCST*F_TxCST-F_TKCST-(F_CSTCNPDL1+F_CSTCPDL1)-F_DCST;

% Non-PDL1, T-exposed tumor cells

F_ResCNPDL1=Res*K_RCNPD*CNPDL1*(1-CNPDL1/Y_CNPDM)*1/(alpha_C*(CST+CPDL1)+1);                          % Resource-driven proliferation of tumor stem cells

F_TxCST=(1+K_TXC*TEX/(1+TEX));                                                              % Exhausted T-cells-driven proliferation modulator of tumor stem cells
 
F_TKCNPDL1= K_TKC*CNPDL1*(TKPD1+TKNPD1);                                                    % T-cell driven apoptosis of the Non-PDL1 tumor cell

F_CNPDL1CPDL1= K_CPDNPD*CNPDL1;                                                             % IFNG-induced conversion from PDL1-ve to +ve tumor cells

F_DCNPDL1= K_CNPD*CNPDL1;                                                                   % Death of Non-PDL1, T-exposed tumor cells

dydt_CNPDL1=F_ResCNPDL1*F_TxCST+F_CSTCNPDL1-F_TKCNPDL1-F_CNPDL1CPDL1-F_DCNPDL1;

% PDL1+, T-exposed tumor cells

F_ResCPDL1=Res*K_RCPD*CPDL1*(1-CPDL1/Y_CPDM)*1/(alpha_C*(CST+CNPDL1)+1);                              % Resource-driven proliferation of PDL1 +ve cells

F_TxCST=(1+K_TXC*TEX/(1+TEX));                                                              % Exhausted T-cells-driven proliferation modulator of tumor stem cells

F_TKNPDCPDL1=K_TKC*CPDL1*TKNPD1;                                                            % Killer, PD1-ve T-cell driven death of PDL1+ve tumor cells

F_DCPDL1= K_CPD*CPDL1;                                                                      % Death of PDL1, T-exposed tumor cells   

dydt_CPDL1=F_ResCPDL1*F_TxCST+F_CNPDL1CPDL1+F_CSTCPDL1-F_TKNPDCPDL1-F_DCPDL1;

%% Resource concentration

F_ProRes=K_RIn*(Y_RM-Res);

F_ConsRes=Res*(K_RCST*CST+K_RCNPD*CNPDL1+K_RCPD*CPDL1);

F_DRes= K_ResD*Res;

dydt_Res=F_ProRes-F_ConsRes-F_DRes;

%% T-cell population

% Killer PD1+ T cells

F_ProTKPD1=K_TKPD*TKPD1*(1-TKPD1/Y_TKPD1);                                             % Proliferation of killer, PD1+ve T-cells

F_TKPD1TKNPD1=K_TKPDNPD1*u*TKPD1;                                                      % ICI-based conversion from PD1+ve to -ve killer T cells

F_TKPD1CPDL1=(1+K_CNPDT*CNPDL1/(CNPDL1+1));                                            % Tumor cells-driven activation of T cells

F_TKPD1TEX=K_TKPDTEX*TKPD1*(CPDL1)/(CPDL1+1);                                          % PDL1+ve and M2 driven conversion of PD1+ve to exhausted T cells  

F_DTKPD1= K_TKPDD*TKPD1;                                                                % Death rate of PD1+, killer T cells

dydt_TKPD1 = F_ProTKPD1*F_TKPD1CPDL1-F_TKPD1TKNPD1-F_TKPD1TEX-F_DTKPD1; 

% Killer PD1- T cells
F_ProTKNPD1=K_TKNPD*TKNPD1*(1-TKNPD1/Y_TKNPD1);                                        % Proliferation of killer, PD1+ve T-cells

F_TKPD1CPDL1=(1+K_CNPDT*CNPDL1/(CNPDL1+1));                                            % Tumor cells-driven activation of T cells

F_DTKNPD1= K_TKNPDD*TKNPD1;

%dydt_TKNPD1=F_TKPD1TKNPD1-F_DTKNPD1;

dydt_TKNPD1=F_ProTKNPD1*F_TKPD1CPDL1+F_TKPD1TKNPD1-F_DTKNPD1;


% Exhausted T cells
F_DTEX=K_TEXD*TEX;                                                                      % Death rate of exhausted t cells

dydt_TEX=F_TKPD1TEX-F_DTEX;




dydt=[dydt_CST;dydt_CNPDL1;dydt_CPDL1;dydt_Res;dydt_TKPD1;dydt_TKNPD1;dydt_TEX];



end

