function [J_AFEiR, J_AFEiX, J_AFErR, J_AFErX] = Jacobian_Converters_1ph(E_re,E_im, E_star, Grid_para, Filter_para, J_EdiorR, J_EdiorX, J_EdioiR, J_EdioiX, J_PR, J_PX, J_QR, J_QX, J_AFEiR, J_AFEiX, J_AFErR, J_AFErX)

    n_dc = Grid_para.n_dc;
    n_ac = Grid_para.n_ac;
    n_nodes = Grid_para.n_nodes;
    n_ph = Grid_para.n_ph;
    n_AFE = Grid_para.n_AFE;
    AFE_type = Grid_para.AFE_type;
    pos_ac3 = Grid_para.pos_ac3;
    pos_dc3 = Grid_para.pos_dc3;
    R = Filter_para.R;
    X = Filter_para.X;
    G = Grid_para.G;
    B = Grid_para.B;
    Y = complex(G,B);
   
   
    %IGBT loss model (accounting for this only improves the accuracy slightly)
    I_b =Grid_para.Y_b*Grid_para.V_b;
       

    for p = 1:n_AFE
        

        if AFE_type(p) == "VdcQ"

        %% Vdc
        
            Imag = abs(Y(pos_ac3(p,1),:)*complex(E_re,E_im));
            R_eq_ctu = interp1(Filter_para.IGBT_piecewise(:,1),Filter_para.IGBT_piecewise(:,2),Imag*I_b)./Grid_para.V_b./Imag;
            R_eq = (R_eq_ctu*4/pi); % ???
            R_eq(isnan(R_eq))=0;
            R_eq(isinf(R_eq))=0;
            
            if Filter_para.Exclude_losses
                R = 0;
                R_eq = 0;
                X = 0;
            end    
            
            
            beta = ((R + R_eq).^2 - X^2) .* G(pos_ac3(p,1),pos_ac3(p,2)) - 2*(R + R_eq).*X.*B(pos_ac3(p,1),pos_ac3(p,2));
            gamma = ((R + R_eq).^2 - X^2) .*B(pos_ac3(p,1),pos_ac3(p,2)) + 2*(R + R_eq).*X.*G(pos_ac3(p,1),pos_ac3(p,2));
            beta = -beta;
            gamma = -gamma;

            P_0pn = (E_re(pos_ac3(p,1)).*( G(pos_ac3(p,1),pos_ac3(p,2))*E_re(pos_ac3(p,2)) + G(pos_ac3(p,1),pos_ac3(p,1))*E_re(pos_ac3(p,1)) - B(pos_ac3(p,1),pos_ac3(p,2))*E_im(pos_ac3(p,2)) - B(pos_ac3(p,1),pos_ac3(p,1))*E_im(pos_ac3(p,1))) + ...
                    E_im(pos_ac3(p,1)).*(B(pos_ac3(p,1),pos_ac3(p,2))*E_re(pos_ac3(p,2)) + B(pos_ac3(p,1),pos_ac3(p,1))*E_re(pos_ac3(p,1)) + G(pos_ac3(p,1),pos_ac3(p,2))*E_im(pos_ac3(p,2)) + G(pos_ac3(p,1),pos_ac3(p,1))*E_im(pos_ac3(p,1)))) + ...
                   (-beta*E_re(pos_ac3(p,2)) + beta*E_re(pos_ac3(p,1)) + gamma*E_im(pos_ac3(p,2)) - gamma*E_im(pos_ac3(p,1)));

               
%             S_loss = (complex(R,X)*A*YY(pos_ac3(p3,1),:)*E_star).*conj(A*YY(pos_ac3(p3,1),:)*E_star)
%             S_ac = A*E_star(pos_ac3(p3,1),:).*conj(A*YY(pos_ac3(p3,1),:)*E_star)   
%             S_dc = E_star(pos_dc3(p,1),:).*conj(YY(pos_dc3(p,1),:)*E_star) 
%                
%             
%             real(( A*E_star(pos_ac3(p3,1),:) + (complex(R,X)*A*YY(pos_ac3(p3,1),:)*E_star)) .*conj(A*YY(pos_ac3(p3,1),:)*E_star)) + S_dc*[0;1;0]
%             real(( A*E_star(pos_ac3(p3,1),:) - (complex(R,X)*A*YY(pos_ac3(p3,1),:)*E_star)) .*conj(A*YY(pos_ac3(p3,1),:)*E_star)) + S_dc*[0;1;0]
%             real(( A*E_star(pos_ac3(p3,1),:) + 0*(complex(R,X)*A*YY(pos_ac3(p3,1),:)*E_star)) .*conj(A*YY(pos_ac3(p3,1),:)*E_star)) + S_dc*[0;1;0]

%             Sdio_star = (Atot*E_star).*conj(Y*(Atot*E_star));
%             P_0pn = real(Sdio_star(pos_ac3(p3,1)))
    
            alpha =  ( G(pos_dc3(p,1),pos_dc3(p,2))* E_re(pos_dc3(p,2)))^2 - 4 * G(pos_dc3(p,1),pos_dc3(p,1)) .* (P_0pn);

            for o = 1:length(alpha)
                if alpha(o) < 1e-12
                    P_0pn
                    alpha
                    alpha(o) = 1
                    disp('Warning: alpha is wrong')
                end
            end

%           Ssysm =   Atot*E_star.*conj(Atot*(G+1i*B)*E_star)
            
            A1_re = -1./sqrt(alpha).* ( G(pos_ac3(p,1) ,pos_ac3(p,2) )*E_re(pos_ac3(p,1) ) + B(pos_ac3(p,1) ,pos_ac3(p,2) )*E_im(pos_ac3(p,1) ) - diag(beta));
            A1_im = -1./sqrt(alpha).* (-B(pos_ac3(p,1) ,pos_ac3(p,2) )*E_re(pos_ac3(p,1) ) + G(pos_ac3(p,1) ,pos_ac3(p,2) )*E_im(pos_ac3(p,1) ) + diag(gamma));

            A2_re = -1./sqrt(alpha).* ( G(pos_ac3(p,1) ,pos_ac3(p,2) )*E_re(pos_ac3(p,2) ) + 2*G(pos_ac3(p,1) ,pos_ac3(p,1) )*E_re(pos_ac3(p,1) ) - B(pos_ac3(p,1) ,pos_ac3(p,2) )*E_im(pos_ac3(p,2) ) + diag(beta));
            A2_im = -1./sqrt(alpha).* ( B(pos_ac3(p,1) ,pos_ac3(p,2) )*E_re(pos_ac3(p,2) ) + G(pos_ac3(p,1) ,pos_ac3(p,2) )*E_im(pos_ac3(p,2) ) + 2*G(pos_ac3(p,1) ,pos_ac3(p,1) )*E_im(pos_ac3(p,1) ) - diag(gamma));


            J_AFErR(pos_ac3(p,1),pos_ac3(p,2) ) = A1_re;
            J_AFErR(pos_ac3(p,1),pos_ac3(p,1) ) = A2_re;
            J_AFErR(pos_ac3(p,1),pos_dc3(p,2) ) = -G(pos_dc3(p,1) ,pos_dc3(p,2) )/(2*G(pos_dc3(p,1) ,pos_dc3(p,1) )) + G(pos_dc3(p,1) ,pos_dc3(p,2) )^2 * E_re(pos_dc3(p,2) ) / ( 2* G(pos_dc3(p,1) ,pos_dc3(p,1) ) * sqrt(alpha) );

            J_AFErX(pos_ac3(p,1),pos_ac3(p,2) ) = A1_im;
            J_AFErX(pos_ac3(p,1),pos_ac3(p,1) ) = A2_im;

        %% Q

        
            J_AFEiR(pos_ac3(p,1),: ) = J_QR(pos_ac3(p,1),: );
            J_AFEiX(pos_ac3(p,1),: ) = J_QX(pos_ac3(p,1),: );

          
        elseif AFE_type(p) == "PQ"

        %% P  

            J_AFErR(pos_ac3(p,1),: ) = J_PR(pos_ac3(p,1),: );
            J_AFErX(pos_ac3(p,1),: ) = J_PX(pos_ac3(p,1),: );
            
            J_AFErR(pos_ac3(p,1),pos_ac3(p,1) ) = J_AFErR(pos_ac3(p,1),pos_ac3(p,1) ) + 2*R*E_re(pos_ac3(p,1))' - 2*R*E_re(pos_ac3(p,2))' - 2*X*E_im(pos_ac3(p,1))' + 2*X*E_im(pos_ac3(p,2))';
            J_AFErR(pos_ac3(p,1),pos_ac3(p,2) ) = J_AFErR(pos_ac3(p,1),pos_ac3(p,2) ) - 2*R*E_re(pos_ac3(p,1))' + 2*R*E_re(pos_ac3(p,2))' + 2*X*E_im(pos_ac3(p,1))' - 2*X*E_im(pos_ac3(p,2))';
            J_AFErX(pos_ac3(p,1),pos_ac3(p,1) ) = J_AFErX(pos_ac3(p,1),pos_ac3(p,1) ) - 2*X*E_re(pos_ac3(p,1))' + 2*X*E_re(pos_ac3(p,2))' - 2*R*E_im(pos_ac3(p,1))' + 2*R*E_im(pos_ac3(p,2))';
            J_AFErX(pos_ac3(p,1),pos_ac3(p,2) ) = J_AFErX(pos_ac3(p,1),pos_ac3(p,2) ) + 2*X*E_re(pos_ac3(p,1))' - 2*X*E_re(pos_ac3(p,2))' + 2*R*E_im(pos_ac3(p,1))' - 2*R*E_im(pos_ac3(p,2))';
           

           
            
            
        %% Q

            J_AFEiR(pos_ac3(p,1),: ) = J_QR(pos_ac3(p,1),: );
            J_AFEiX(pos_ac3(p,1),: ) = J_QX(pos_ac3(p,1),: );

      
        end
    end

end