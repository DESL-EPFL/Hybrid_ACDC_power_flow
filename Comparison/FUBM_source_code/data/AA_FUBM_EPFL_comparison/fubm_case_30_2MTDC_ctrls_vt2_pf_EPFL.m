function mpc = fubm_case_30_2MTDC_ctrls_vt2_pf_EPFL
%%%%%%%%%%%%%%%%%%%%%%%% IEEE_30_VSC  Test Case %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% area data
% area refbus
mpc.areas = [
        1       101; %area 1 IEEE 30 Bus System
        2      1001; %area 2 DC GRID 1
        3      4001; %area 3 DC GRID 1
                    %area 3 Transformers for VSC
            ];
        

%% bus data
%	bus_i	type	Pd      Qd      Gs      Bs      area   Vm       Va    base  	zone	Vmax	Vmin

mpc.bus = [
    %Area 1 (IEEE 30 Bus System)
    101     3       0       0       0       0       1       1       0       135     1       1.15	0.9;
    102     2   	30  	12.7	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    103 	1   	2.4     1.2 	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    104 	1   	7.6 	1.6     0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    105 	1   	0   	0   	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    106 	1   	0   	0   	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    107 	1   	22.8	10.9	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    108 	1   	30  	30  	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    109 	1   	0   	0   	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    110 	1   	5.8 	2   	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    111 	1   	0   	0   	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    112 	1   	11.2	7.5 	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    113 	2   	-37   	-15   	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    114 	1   	6.2 	1.6 	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    115 	1   	8.2 	2.5 	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    116 	1   	3.5 	1.8 	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    117 	1   	9   	5.8 	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    118 	1   	3.2 	0.9 	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    119 	1   	9.5 	3.4 	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    120 	1   	2.2 	0.7 	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    121 	1   	17.5	11.2	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    122 	2   	0   	0   	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    123 	2   	3.2 	1.6 	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    124 	1   	8.7 	6.7 	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    125 	1   	0   	0   	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    126 	1   	3.5 	2.3 	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    127 	2   	0   	0   	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    128 	1   	0   	0   	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    129 	1   	2.4 	0.9 	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    130 	1   	10.6	1.9 	0   	0   	1   	1   	0   	135 	1   	1.15	0.9;
    %Area 2 (Transformers nodes for VSC DC GRID 1)  
    901 	1   	0   	0   	0   	0		4   	1   	0   	200 	1   	1.15	0.9;
    902 	1   	0   	0   	0   	0		4   	1   	0   	200 	1   	1.15	0.9;
    903 	1   	0   	0   	0   	0		4   	1   	0   	200 	1   	1.15	0.9;
    %Area 3 (Transformers nodes for VSC DC GRID 1)     
    904 	1   	0   	0   	0   	0		5   	1   	0   	200 	1   	1.15	0.9;
    905 	1   	0   	0   	0   	0		5   	1   	0   	200 	1   	1.15	0.9;
    906 	1   	0   	0   	0   	0		5   	1   	0   	200 	1   	1.15	0.9
    %907 	1   	0   	0   	0   	0		5   	1   	0   	200 	1   	1.15	0.9;
   % 908 	1   	0   	0   	0   	0		5   	1   	0   	200 	1   	1.15	0.9
    %Area 2 (DC GRID 1)
    1    	1   	0   	0   	0   	0   	2   	1   	0   	200 	1   	1.15	0.9;
    2   	1   	0   	0   	0   	0   	2   	1   	0   	200 	1   	1.15	0.9;
    3   	1   	0   	0   	0   	0   	2   	1   	0   	200 	1   	1.15	0.9;
    %Area 3 (DC GRID 2)
    4   	1   	0   	0   	0   	0   	3   	1   	0   	200 	1   	1.15	0.9;
    5   	1   	0   	0   	0   	0   	3   	1   	0   	200 	1   	1.15	0.9;
    6   	1   	0   	0   	0   	0   	3   	1   	0   	200 	1   	1.15	0.9;
    7   	1   	-15 	0   	0   	0   	3   	1   	0   	200 	1   	1.15	0.9;
    8   	1   	-10 	0   	0   	0   	3   	1   	0   	200 	1   	1.15	0.9;
    %Area 2 extra (DC GRID 1)
    1001    1   	0   	0   	0   	0   	6   	1   	0   	200 	1   	1.15	0.9;
    2001    1   	0   	0   	0   	0   	6   	1   	0   	200 	1   	1.15	0.9;
    3001    1   	0   	0   	0   	0   	6   	1   	0   	200 	1   	1.15	0.9;
    %Area 3 extra (DC GRID 2)
    4001    1   	0   	0   	0   	0   	7   	1   	0   	200 	1   	1.15	0.9;
    5001    1   	0   	0   	0   	0   	7   	1   	0   	200 	1   	1.15	0.9;
    6001    1   	0   	0   	0   	0   	7   	1   	0   	200 	1   	1.15	0.9;
   % 7001    1   	0   	0   	0   	0   	7   	1   	0   	200 	1   	1.15	0.9;
   % 8001    1   	0   	0   	0   	0   	7   	1   	0   	200 	1   	1.15	0.9;
    

];
%% generator data
%	bus     Pg      Qg      Qmax	Qmin	Vg      mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
    %IEEE 57 Bus System
    101 	23.54	0   	150 	-20 	1   	100 	1   	80  	0   	0	0	0   	0   	0   	0   	0       	0   	0   	0   	0;
    102 	60.97	0   	60  	-20 	1   	100 	1   	80  	0   	0	0	0   	0   	0   	0   	0       	0   	0   	0   	0;
    122 	21.59	0   	62.5	-15 	1   	100 	1   	50  	0   	0	0	0   	0   	0   	0   	0       	0   	0   	0   	0;
    127 	26.91	0   	48.7	-15 	1   	100 	1   	55  	0   	0	0	0   	0   	0   	0   	0       	0   	0   	0   	0;
    123 	19.2	0   	40      -10 	1   	100 	1   	30  	0   	0	0	0   	0   	0   	0   	0       	0   	0   	0   	0;
%     113 	37  	0   	44.7    -15 	1   	100 	1   	40  	0   	0	0	0   	0   	0   	0   	0       	0   	0   	0   	0
    %WF Generator is added as a load in bus 7 and 8
    ];
%% branch data
%  fbus tbus	r	    x	    b           rateA   rateB	rateC ratio/ma  angle status angmin	angmax  PF      QF      PT      QT      MU_SF  MU_ST   MU_ANGMIN   MU_ANGMAX   VF_SET  VT_SET  MA_MAX MA_MIN CONV_A     BEQ     K2   BEQ_MIN BEQ_MAX SH_MIN  SH_MAX  GSW ALPHA1 ALPHA2 ALPHA3  KDP ----------------------
mpc.branch = [
%Area 1  IEEE 57 Bus System  
    101	102 	0.02	0.06	0.03        130 	130 	130 	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0       1       1       0	    0       1    0       0       -360    360     0    0       0       0     0;
    101	103 	0.05	0.19	0.02        130 	130 	130 	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0         	0   	0       1       1       0   	0   	1    0       0       -360    360     0    0       0       0     0;
    102	104 	0.06	0.17	0.02    	65  	65  	65  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0         	0   	0   	1   	1       0   	0   	1    0       0       -360    360     0    0       0       0     0;
    103	104 	0.01	0.04	0           130 	130 	130 	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0          	0           0   	0   	1   	1       0   	0   	1    0       0       -360    360     0    0       0       0     0;
    102	105 	0.05	0.2 	0.02        130 	130 	130 	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    102	106 	0.06	0.18	0.02    	65      65  	65  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    104	106 	0.01	0.04	0       	90  	90  	90  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    105	107 	0.05	0.12	0.01    	70  	70  	70  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    106	107 	0.03	0.08	0.01    	130 	130 	130 	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    106	108 	0.01	0.04	0           32  	32  	32  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    106	109 	0   	0.21	0           65  	65  	65  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    106	110 	0   	0.56	0           32  	32  	32  	1   	0   	1	-360	360 	0   	0   	0   	3.0   	0   	0   	0       	0       	0   	0   	1.2   	0.8   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    109	111 	0   	0.21	0           65  	65  	65  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    109	110 	0   	0.11	0           65  	65  	65  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    104	112 	0   	0.26	0           65  	65  	65  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    112	113 	0   	0.14	0          	65  	65  	65  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    112	114 	0.12	0.26	0       	32  	32  	32  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0      -50.0   50.0     0    0       0       0     0;
    112	115 	0.07	0.13	0       	32  	32  	32  	1   	0   	1	-360	360 	0     	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    112	116 	0.09	0.2 	0       	32  	32  	32  	1   	0   	1	-360	360 	0       0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    114	115 	0.22	0.2 	0       	16  	16  	16  	1   	0   	1	-360	360     0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    116	117 	0.08	0.19	0       	16  	16  	16  	1   	0   	1	-360	360 	0    	0   	0   	0   	0   	0   	0       	0       	0   	1.01 	1.2   	0.8   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    115	118 	0.11	0.22	0       	16  	16  	16  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    118	119 	0.06	0.13	0       	16  	16  	16  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    119	120 	0.03	0.07	0       	32  	32  	32  	1   	0   	1	-360	360 	0   	0   	0   	0       0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    110	120 	0.09	0.21	0       	32  	32  	32  	1   	0   	1	-360	360 	0   	0       0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0      -50.0   50.0     0    0       0       0     0;%PST1
    110	117  	0.03	0.08	0       	32  	32  	32  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    110	121 	0.03	0.07	0       	32  	32  	32  	1   	0       1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    110	122 	0.07	0.15	0       	32  	32  	32  	1   	0   	1	-360	360 	0       0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    121	122 	0.01	0.02	0       	32  	32  	32  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    115	123 	0.1 	0.2 	0       	16  	16  	16  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    122	124 	0.12	0.18	0       	16  	16  	16  	1   	0   	1	-360	360 	0   	0   	0   	0   	0       0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    123	124 	0.13	0.27	0       	16  	16  	16  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    124	125 	0.19	0.33	0       	16  	16  	16  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0       0   	1    0       0       -360    360     0    0       0       0     0;
    125	126 	0.25	0.38	0       	16  	16  	16  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    125	127 	0.11	0.21	0       	16  	16  	16  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    128	127 	0   	0.4 	0       	65  	65  	65  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0       1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    127	129 	0.22	0.42	0       	16  	16  	16  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    127	130 	0.32	0.6 	0       	16  	16  	16  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    129	130 	0.24	0.45	0       	16  	16  	16  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    108	128 	0.06	0.2 	0.02    	32  	32  	32  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0           0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    106	128 	0.02	0.06	0.01    	32  	32  	32  	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
%Area 1-2  (Transformers for VSC DC GRID 1) rateA   rateB	rateC ratio/ma  angle status angmin	angmax  PF      QF      PT      QT      MU_SF  MU_ST   MU_ANGMIN   MU_ANGMAX   VF_SET  VT_SET  MA_MAX MA_MIN CONV_A     BEQ     K2   BEQ_MIN BEQ_MAX SH_MIN  SH_MAX  GSW ALPHA1 ALPHA2 ALPHA3  ----------------------
    113	901 	0.0015	0.1121	0       	100 	100 	100 	1   	0    	1	-360	360 	0   	0   	0   	0   	0   	0   	0           0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    115	902 	0.0015	0.1121	0       	100 	100 	100 	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0           0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    130	903 	0.0015	0.1121	0       	100 	100 	100 	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0           0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
%Area 1-3  (Transformers for VSC DC GRID 2) rateA   rateB	rateC ratio/ma  angle status angmin	angmax  PF      QF      PT      QT      MU_SF  MU_ST   MU_ANGMIN   MU_ANGMAX   VF_SET  VT_SET  MA_MAX MA_MIN CONV_A     BEQ     K2   BEQ_MIN BEQ_MAX SH_MIN  SH_MAX  GSW ALPHA1 ALPHA2 ALPHA3  ----------------------     
    102	904 	0.0015	0.1121	0       	100 	100 	100 	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    105	905 	0.0015	0.1121	0       	100 	100 	100 	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    106	906 	0.0015	0.1121	0       	100 	100 	100 	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0           0       	0   	0       1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
%Area 2  (VSC of DC GRID 1)                 rateA   rateB	rateC ratio/ma  angle status angmin	angmax  PF      QF      PT      QT      MU_SF  MU_ST   MU_ANGMIN   MU_ANGMAX   VF_SET  VT_SET  MA_MAX MA_MIN CONV_A     BEQ     K2   BEQ_MIN BEQ_MAX SH_MIN  SH_MAX  GSW ALPHA1 ALPHA2 ALPHA3  ----------------------    
   1001 901 	0.0001	0.1643	0       	100 	100 	100 	1   	0   	1	-360	360 	0   	0   	0   	-0.2   	0   	0       0       	0       	1   	0   	1.4     0.6    	2   	0    	1    -0.5    0.5     -360    360     0    0       0       0     0; %VSC1vdc
   2001	902 	0.0001	0.1643	0       	100 	100 	100 	1   	0   	1	-360	360 	0   	0   	0   	-0.2   	0   	0   	0       	0       	1   	0   	1.4     0.6     2   	0	 	1    -0.5    0.5     -50     50      0    0       0       0     0; %VSC2z 
   3001 903 	0.0001	0.1643	0       	100 	100 	100 	1   	0   	1	-360	360 	0.1  	0   	0   	-0.2   	0   	0   	0       	0       	0   	0   	1.4     0.6    	1   	0	 	1    -0.5    0.5     -50     50      0    0       0       0     0; %VSC3z 
%Area 3  (VSC of DC GRID 2)                 rateA   rateB	rateC ratio/ma  angle status angmin	angmax  PF      QF      PT      QT      MU_SF  MU_ST   MU_ANGMIN   MU_ANGMAX   VF_SET  VT_SET  MA_MAX MA_MIN CONV_A     BEQ     K2   BEQ_MIN BEQ_MAX SH_MIN  SH_MAX  GSW ALPHA1 ALPHA2 ALPHA3  ----------------------     
   4001	904 	0.0001	0.1643	0       	100 	100 	100 	1   	0   	1	-360	360 	0   	0   	0   	-0.2   	0   	0   	0           0       	1   	0   	1.4     0.6    	2   	0	    1    -0.5    0.5     -360    360     0    0       0       0     0; %VSC4vdc
   5001	905 	0.0001	0.1643	0       	100 	100 	100 	1       0   	1	-360	360 	0.1    	0   	0   	-0.2   	0   	0   	0       	0           0   	0       1.4     0.6    	1   	0	 	1    -0.5    0.5     -50     50      0    0       0       0     0; %VSC5z 
   6001	906 	0.0001	0.1643	0       	100 	100 	100 	1   	0   	1	-360	360 	0.1   	0   	0   	-0.2   	0   	0   	0           0       	0   	0   	1.4     0.6    	1   	0	 	1    -0.5    0.5     -50     50      0    0       0       0     0; %VSC6z 
   %  7001	907 	0.0001	0.1643	0       	100 	100 	100 	1       0   	1	-360	360 	0.1    	0   	0   	0.2   	0   	0   	0       	0           0   	0       1.4     0.6    	1   	0	 	1    -0.5    0.5     -50     50      0    0       0       0     0; %VSC7z 
 %  8001	908 	0.0001	0.1643	0       	100 	100 	100 	1   	0   	1	-360	360 	0.1   	0   	0   	0.2   	0   	0   	0           0       	0   	0   	1.4     0.6    	1   	0	 	1    -0.5    0.5     -50     50      0    0       0       0     0; %VSC8vdc 
    %Area 2  (DC GRID 1)                        rateA   rateB	rateC ratio/ma  angle status angmin	angmax  PF      QF      PT      QT      MU_SF  MU_ST   MU_ANGMIN   MU_ANGMAX   VF_SET  VT_SET  MA_MAX MA_MIN CONV_A     BEQ     K2   BEQ_MIN BEQ_MAX SH_MIN  SH_MAX  GSW ALPHA1 ALPHA2 ALPHA3  ----------------------   
    001	002 	0.05	0   	0        	200 	200 	200     1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0           0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    001	003 	0.05	0   	0       	200 	200 	200 	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    002	003 	0.05	0   	0       	200 	200 	200 	1   	0   	1	-360	360 	0   	0      	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    001 1001    1e-3    0       0           200     200     200     1       0       1   -360	360     0   	0      	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    002 2001    1e-3    0       0           200     200     200     1       0       1   -360	360     0   	0      	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    003 3001    1e-3    0       0           200     200     200     1       0       1   -360	360     0   	0      	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
%Area 3  (DC GRID 2)                        rateA   rateB	rateC ratio/ma  angle status angmin	angmax  PF      QF      PT      QT      MU_SF  MU_ST   MU_ANGMIN   MU_ANGMAX   VF_SET  VT_SET  MA_MAX MA_MIN CONV_A     BEQ     K2   BEQ_MIN BEQ_MAX SH_MIN  SH_MAX  GSW ALPHA1 ALPHA2 ALPHA3  ----------------------     
    004	005     0.05	0   	0       	200 	200 	200 	1   	0   	1	-360	360 	0   	0   	0   	0    	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    004	006 	0.05	0   	0       	200 	200 	200 	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    005	008 	0.05	0   	0       	200 	200 	200 	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    006	007 	0.05	0   	0       	200 	200 	200 	1   	0   	1	-360	360 	0   	0   	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    004 4001    1e-3    0       0           200     200     200     1       0       1   -360	360     0   	0      	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    005 5001    1e-3    0       0           200     200     200     1       0       1   -360	360     0   	0      	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
    006 6001    1e-3    0       0           200     200     200     1       0       1   -360	360     0   	0      	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
 %   007 7001    1e-5    0       0           200     200     200     1       0       1   -360	360     0   	0      	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
 %   008 8001    1e-5    0       0           200     200     200     1       0       1   -360	360     0   	0      	0   	0   	0   	0   	0       	0       	0   	0   	1   	1   	0   	0   	1    0       0       -360    360     0    0       0       0     0;
];
%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
    %Pg cost
    2	0   	0       	3	0.02        	2   	0;
    2	0   	0       	3	0.0175      	1.75	0;
    2	0   	0       	3	0.0625      	1   	0;
    2	0   	0       	3	0.00834     	3.25	0;
    2	0   	0       	3	0.025       	3   	0;
    2	0   	0       	3	0.025       	3   	0;
    %Qg cost
    2	0   	0       	3	0.0000        	0   	0;
    2	0   	0       	3	0.0000        	0   	0;
    2	0   	0       	3	0.0000        	0   	0;
    2	0   	0       	3	0.0000        	0   	0;
    2	0   	0       	3	0.0000        	0   	0;
    2	0   	0       	3	0.0000        	0   	0;

];
