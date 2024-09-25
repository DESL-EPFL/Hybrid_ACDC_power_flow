function mpc = FUBM_test_grid
%%%%%%%%%%%%%%%%%%%%%%%% IEEE_30_VSC  Test Case %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 0.1;

%% area data
% area refbus
mpc.areas = [
        1    1; %area 1 IEEE 30 Bus System
        2    4; %area 2 DC GRID 1
            ];

        

%% bus data
%	bus_i	type	Pd      Qd      Gs      Bs      area   Vm       Va    base[KV]	zone	Vmax	Vmin

mpc.bus = [
    %Area 1 (IEEE 30 Bus System)
    1       3       0       0       0       0       1       1       0       0.4     1       1.05	0.95;
    2       1   	0.0075 -0.015   0   	0   	1   	1   	0   	0.4 	1   	1.05	0.95;
    3    	1   	0       0   	0   	0    	1   	1   	0   	0.4 	1   	1.05	0.95;
    %Area 2 (Transformers nodes for VSC DC GRID 1)  
    4    	1   	0   	0   	0   	0   	2   	1   	0   	0.8 	2   	1.05	0.95;
    %Area 2 (DC GRID 1)
    5   	1   	0.01    0   	0   	0   	2   	1   	0   	0.8 	2   	1.05	0.95;
];
%% generator data
%	bus     Pg      Qg      Qmax	Qmin	Vg     mBase	status	Pmax   Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
    %IEEE 57 Bus System
    1       0       0   	1       -1  	1   	0.1    	1   	1  	   -1   0	0	0   	0   	0   	0   	0       	0   	0   	0   	0;
    ];
%% branch data
%  fbus tbus	r	    x	    b       rateA   rateB	rateC ratio/ma  angle status angmin	angmax  PF      QF      PT      QT      MU_SF  MU_ST   MU_ANGMIN   MU_ANGMAX   VF_SET  VT_SET  MA_MAX MA_MIN CONV_A     BEQ     K2   BEQ_MIN BEQ_MAX SH_MIN  SH_MAX  GSW ALPHA1 ALPHA2 ALPHA3  KDP----------------------
mpc.branch = [
%Area 1  IEEE 57 Bus System  
    1	2 	0.0051	0.0022	0           0.1 	0.1 	0.1 	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0       1       1       0	    0       1    0       0       -360    360     0    0       0       0     0;
%Area 1-2  (Transformers for VSC DC GRID 1) rateA   rateB	rateC ratio/ma  angle status angmin	angmax  PF      QF      PT      QT      MU_SF  MU_ST   MU_ANGMIN   MU_ANGMAX   VF_SET  VT_SET  MA_MAX MA_MIN CONV_A     BEQ     K2   BEQ_MIN BEQ_MAX SH_MIN  SH_MAX  GSW ALPHA1 ALPHA2 ALPHA3  ----------------------
    2	3 	0.0051	0.0022	0       	0.1 	0.1 	0.1 	1   	0    	1	-360	360 	0   	0   	0   	0   	0   	0   	0           0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
%Area 2  (VSC of DC GRID 1)                 rateA   rateB	rateC ratio/ma  angle status angmin	angmax  PF      QF      PT      QT      MU_SF  MU_ST   MU_ANGMIN   MU_ANGMAX   VF_SET  VT_SET  MA_MAX MA_MIN CONV_A     BEQ     K2   BEQ_MIN BEQ_MAX SH_MIN  SH_MAX  GSW ALPHA1 ALPHA2 ALPHA3  ----------------------    
    4	3 	0.0050  0.025   0       	0.1 	0.1 	0.1 	1   	0   	1	-360	360 	0   	0   	0   -0.01   	0   	0       -1       	1       	1   	0       1.2     0.8     2       0       1    -0.7    0.7     -360    360     0    0       0       0     0; %VSC1vdc
%Area 2  (DC GRID 1)                        rateA   rateB	rateC ratio/ma  angle status angmin	angmax  PF      QF      PT      QT      MU_SF  MU_ST   MU_ANGMIN   MU_ANGMAX   VF_SET  VT_SET  MA_MAX MA_MIN CONV_A     BEQ     K2   BEQ_MIN BEQ_MAX SH_MIN  SH_MAX  GSW ALPHA1 ALPHA2 ALPHA3  ----------------------   
    4	5 	0.0625	0   	0        	0.1 	0.1 	0.1     1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0           0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
];
%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
    %Pg cost
    2	0   	0       	1	0.01        	2   	0;
    %Qg cost
    2	0   	0       	1	0.01        	0   	0;
];

