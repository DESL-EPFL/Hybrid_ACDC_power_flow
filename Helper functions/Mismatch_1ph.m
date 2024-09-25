function dF = Mismatch_1ph(E,S,E_star,S_star,Grid_para,Filter_para,idx)
    
    n_dc = Grid_para.n_dc;
    n_ac = Grid_para.n_ac;
    n_nodes = Grid_para.n_nodes;
    n_ph = Grid_para.n_ph;
    G = Grid_para.G;
    B = Grid_para.B;
    Y = complex(G,B);    
    R = Filter_para.R;
    X = Filter_para.X;
    I_b =Grid_para.Y_b*Grid_para.V_b;
    
    pos_ac3 = Grid_para.pos_ac3;
    pos_dc3 = Grid_para.pos_dc3;
    
    % Recompute the DC voltage at the AFE
    
    for p = 1:length(idx.vscdc_vq)
        alphan =  ( G(pos_dc3(p,1),pos_dc3(p,2))* E(pos_dc3(p,2)) )^2 - 4 * G(pos_dc3(p,1),pos_dc3(p,1)) * real( S(idx.vscac_vq(p)));
        E(idx.vscdc_vq(p)) = -G(pos_dc3(p,1),pos_dc3(p,2))/(2*G(pos_dc3(p,1),pos_dc3(p,1))).* E(pos_dc3(p,2)) + sqrt(alphan)/(2*G(pos_dc3(p,1),pos_dc3(p,1)));
    end
        
    Imag = abs(Y*E);
    R_eq_ctu = interp1(Filter_para.IGBT_piecewise(:,1),Filter_para.IGBT_piecewise(:,2),Imag*I_b)./Grid_para.V_b./Imag;
    R_eq = (R_eq_ctu*4/pi); % ???
    R_eq(isnan(R_eq))=0;
    R_eq(isinf(R_eq))=0;
            
    if Filter_para.Exclude_losses
        R = 0;
        R_eq = 0;
        X = 0;
    end    

    % include the voltage drop over the filter AND IGBT
    E_filter = E + ((R + R_eq) + 1i*X) .* complex(G,B) * E;
    E(idx.vscdc_vq) = E_filter(idx.vscdc_vq);

    
    % Compute the mismatches for the entire network.
    dS = S_star-S;
    dP = real(dS);
    dQ = imag(dS);
    dV2 = abs(E_star).^2-abs(E).^2; % assumption: V^{*}=1 
    dPdc = real(dS);
    dEdc = real(E_star-E);

    
    
    dF = [ dP(idx.pqac); 
           dQ(idx.pqac); 
           dP(idx.pvac); 
           dV2(idx.pvac);              
           dP(idx.vscac_pq);
           dQ(idx.vscac_pq);
           dEdc(idx.vscdc_vq);
           dQ(idx.vscac_vq);
           dPdc(idx.vscdc_pq);
           dPdc(idx.pdc)];
     
end