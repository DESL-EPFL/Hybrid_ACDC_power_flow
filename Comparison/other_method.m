A_b = 1e5;
V_b= 400;
Y_b = A_b/V_b^2; 
I_b=A_b/(V_b*sqrt(3));

Vdc_b = 800;
Adc_b = A_b;
Ydc_b = Adc_b/Vdc_b^2; 
Idc_b = Adc_b/Vdc_b;

Filter_para.R = 0.016*Y_b; %checked 
Filter_para.X = 0.000127325*2*pi*50*Y_b;  %checked


Filter_para.IGBT_piecewise = [  0                   0
                                0.04926559627563	0.7
                                2.30625399327864	0.8
                                15.7793399043317	0.85
                                107.547461516782	0.8999
                                735.837403888342	0.9499
                                1588.01477341768	0.9699];
Filter_para.Include_losses = 1;


data = load('/Users/willem/Documents/phd/State_estimation/EMTP/Experiments/data_SC_4.mat'); % balanced + filter


%% VSI 1
V1 = transpose(complex(data.B15_Va_mag_control.*cos(data.B15_Va_rad_control), data.B15_Va_mag_control.*sin(data.B15_Va_rad_control)))/(sqrt(2)/sqrt(3)*V_b);
V1ppp = transpose(complex(data.VSI1_Bppp_Va_mag_control.*cos(data.VSI1_Bppp_Va_rad_control), data.VSI1_Bppp_Va_mag_control.*sin(data.VSI1_Bppp_Va_rad_control)))/(sqrt(2)/sqrt(3)*V_b);
V1pp = transpose(complex(data.VSI1_Bpp_Va_mag_control.*cos(data.VSI1_Bpp_Va_rad_control), data.VSI1_Bpp_Va_mag_control.*sin(data.VSI1_Bpp_Va_rad_control)))/(sqrt(2)/sqrt(3)*V_b);
V1p = transpose(complex(data.VSI1_Bp_Va_mag_control.*cos(data.VSI1_Bp_Va_rad_control), data.VSI1_Bp_Va_mag_control.*sin(data.VSI1_Bp_Va_rad_control)))/(sqrt(2)/sqrt(3)*V_b);
I1 = transpose(complex(data.VSI1_Ia_mag_control.*cos(data.VSI1_Ia_rad_control), data.VSI1_Ia_mag_control.*sin(data.VSI1_Ia_rad_control)))/sqrt(2)./I_b;
Vdc1 = transpose(data.B19_Vdc_flow_p_control)/Vdc_b;
Idc1 = transpose(data.B19_Idc_flow_p_control)/Idc_b;
M1 = transpose(complex(data.VSI1_Mreal_a_2_control,data.VSI1_Mimag_a_2_control));
V1 = V1(1000:end);
V1ppp = V1ppp(1000:end);
V1pp = V1pp(1000:end);
V1p = V1p(1000:end);
I1 = I1(1000:end);
Vdc1 = Vdc1(1000:end);
Idc1 = Idc1(1000:end);
M1 = M1(1000:end);

%% VSI 3
V3 = transpose(complex(data.B17_Va_mag_control.*cos(data.B17_Va_rad_control), data.B17_Va_mag_control.*sin(data.B17_Va_rad_control)))/(sqrt(2)/sqrt(3)*V_b);
V3ppp = transpose(complex(data.VSI3_Bppp_Va_mag_control.*cos(data.VSI3_Bppp_Va_rad_control), data.VSI3_Bppp_Va_mag_control.*sin(data.VSI3_Bppp_Va_rad_control)))/(sqrt(2)/sqrt(3)*V_b);
V3pp = transpose(complex(data.VSI3_Bpp_Va_mag_control.*cos(data.VSI3_Bpp_Va_rad_control), data.VSI3_Bpp_Va_mag_control.*sin(data.VSI3_Bpp_Va_rad_control)))/(sqrt(2)/sqrt(3)*V_b);
V3p = transpose(complex(data.VSI3_Bp_Va_mag_control.*cos(data.VSI3_Bp_Va_rad_control), data.VSI3_Bp_Va_mag_control.*sin(data.VSI3_Bp_Va_rad_control)))/(sqrt(2)/sqrt(3)*V_b);
I3 = transpose(complex(data.VSI3_Ia_mag_control.*cos(data.VSI3_Ia_rad_control), data.VSI3_Ia_mag_control.*sin(data.VSI3_Ia_rad_control)))/sqrt(2)./I_b;
Vdc3 = transpose(data.B21_Vdc_flow_p_control)/Vdc_b;
Idc3 = transpose(data.B21_Idc_flow_p_control)/Idc_b;

V3 = V3(1000:end);
V3ppp = V3ppp(1000:end);
V3pp = V3pp(1000:end);
V3p = V3p(1000:end);
I3 = I3(1000:end);
R_eq3 = interp1(Filter_para.IGBT_piecewise(:,1),Filter_para.IGBT_piecewise(:,2),abs(I3)*sqrt(2)*I_b)./V_b./(abs(I3)/sqrt(3)*sqrt(2))*4/pi;


%% 2-port model from literature

Gsw = 0;
ys = 0;
bc = 
Beq = 
ma = abs(M1);
theta = angle(M1);

Y11 = Gsw  + (ys + 1i*bc/2 + 1i*Beq)/ma^2;
Y12 = -ys/(ma*exp(-1i*theta));
Y21 = -ys/(ma*exp(1i*theta));
Y22 = ys + 1i*bc/2;

Y = [Y11 Y12;
     Y21 Y22];



