function dF = Mismatch(E,S,E_star,S_star,Grid_para,Filter_para,idx)
    
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
    
    alp = exp(2*pi/3*1i);
    A = 1/3*[1 1     1; 
             1 alp   alp^2; 
             1 alp^2 alp];

    ACell =  repmat({A}, 1, n_ac);
    ICell =  repmat({1}, 1, n_dc); 
    Atot = blkdiag(ACell{:},ICell{:});
    
    
    % Recompute the DC voltage at the AFE
%     Sdio_star = (Atot*E_star).*conj(Y*(Atot*E_star));
    Sdio = (Atot*E).*conj(Atot*Y*(E));
%     Sdio_delta = Sdio_star - Sdio;
    
    for p = 1:length(idx.vscdc_vq)
        
        p3 = polyphase_indices(p,3);
        alphan =  ( G(pos_dc3(p,1),pos_dc3(p,2))* E(pos_dc3(p,2)) )^2 - 4 * G(pos_dc3(p,1),pos_dc3(p,1)) * real( Sdio(idx.vscac_vq(p3(2))));
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
    dPp =  real((Atot*E_star).*conj(Y*(Atot*E_star)) - (Atot*E).*conj(Y*(Atot*E)));
    dQp =  imag((Atot*E_star).*conj(Y*(Atot*E_star)) - (Atot*E).*conj(Y*(Atot*E)));
    dEdior = real(Atot*E_star) - real(Atot*E_filter);
    dEdioi = imag(Atot*E_star) - imag(Atot*E_filter);
    
    
    dF = [ dP(idx.pqac); 
           dQ(idx.pqac); 
           dP(idx.pvac); 
           dV2(idx.pvac);              
           matrix2vector([ dEdior(idx.vscac_pq(1:3:end))';
                           dPp(idx.vscac_pq(2:3:end))';
                           dEdior(idx.vscac_pq(3:3:end))']);           
           matrix2vector([ dEdioi(idx.vscac_pq(1:3:end))';
                           dQp(idx.vscac_pq(2:3:end))';
                           dEdioi(idx.vscac_pq(3:3:end))']);   
           matrix2vector([ dEdior(idx.vscac_vq(1:3:end))';
                           dEdc(idx.vscdc_vq)';
                           dEdior(idx.vscac_vq(3:3:end))']);           
           matrix2vector([ dEdioi(idx.vscac_vq(1:3:end))';
                           dQp(idx.vscac_vq(2:3:end))';
                           dEdioi(idx.vscac_vq(3:3:end))']);    
           dPdc(idx.vscdc_pq);
           dPdc(idx.pdc)];
           %dEdc(idx.vscdc_vq)];
        
end