function data = process_EMTP_data(emtp_data)

    [Nodal_V_mag,Nodal_V_angle, Nodal_I_mag, Nodal_I_angle, Flow_I_mag, Flow_I_angle, Idc_flow, Idc_inj, Vdc_LF, n_timesteps, Mreal, Mimag, Mabs,Nodal_P,Nodal_Q,Pdc_inj,M_LF] = GetEMTPdata_LF(emtp_data,A_b,V_b,Adc_b, Vdc_b,1,[],Grid_para.n_ph);

    Nodal_V_mag = Nodal_V_mag(end,:);
    Nodal_V_angle =  Nodal_V_angle(end,:);

    Nodal_P =  Nodal_P(end,:);
    Nodal_Q =  Nodal_Q(end,:);
    Pdc_inj =  Pdc_inj(end,:);

    V_complex_LF = transpose(complex(Nodal_V_mag.*cos(Nodal_V_angle), Nodal_V_mag.*sin(Nodal_V_angle)));
    Vdc_LF = transpose(Vdc_LF);

    data.E_star = [V_complex_LF(1:end,1); Vdc_LF(:,1)];
    data.S_star = [transpose(complex(Nodal_P, Nodal_Q)); transpose(complex(Pdc_inj,0))];

end