function [J_AFEiR, J_AFEiX, J_AFErR, J_AFErX] = Jacobian_Converters(E_re,E_im, E_star, Grid_para, Filter_para, J_EdiorR, J_EdiorX, J_EdioiR, J_EdioiX, J_PpR, J_PpX, J_QpR, J_QpX, J_AFEiR, J_AFEiX, J_AFErR, J_AFErX)

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
   
    alp = exp(2*pi/3*1i);
    A = 1/3*[1 1     1; 
             1 alp   alp^2; 
             1 alp^2 alp];

    ACell =  repmat({A}, 1, n_ac);
    ICell =  repmat({1}, 1, n_dc); 
    Atot = blkdiag(ACell{:},ICell{:});

    Edio_re = real(Atot*(complex(E_re,E_im)));
    Edio_im = imag(Atot*(complex(E_re,E_im)));

    %IGBT loss model (accounting for this only improves the accuracy slightly)
    I_b =Grid_para.Y_b*Grid_para.V_b;
       

    for p = 1:n_AFE
        
        p3 = polyphase_indices(p,3);

        if AFE_type(p) == "VdcQ"

        %% Vdc
        
            Imag = abs(Y(pos_ac3(p3,1),:)*complex(E_re,E_im));
            R_eq_ctu = interp1(Filter_para.IGBT_piecewise(:,1),Filter_para.IGBT_piecewise(:,2),Imag*I_b)./Grid_para.V_b./Imag;
            R_eq = (R_eq_ctu*4/pi); % ???
            R_eq(isnan(R_eq))=0;
            R_eq(isinf(R_eq))=0;
            
            if Filter_para.Exclude_losses
                R = 0;
                R_eq = 0;
                X = 0;
            end    
            
            
            beta = ((R + R_eq).^2 - X^2) .* G(pos_ac3(p3,1),pos_ac3(p3,2)) - 2*(R + R_eq).*X.*B(pos_ac3(p3,1),pos_ac3(p3,2));
            gamma = ((R + R_eq).^2 - X^2) .*B(pos_ac3(p3,1),pos_ac3(p3,2)) + 2*(R + R_eq).*X.*G(pos_ac3(p3,1),pos_ac3(p3,2));
            beta = -beta;
            gamma = -gamma;

            P_0pn = (Edio_re(pos_ac3(p3,1)).*( G(pos_ac3(p3,1),pos_ac3(p3,2))*Edio_re(pos_ac3(p3,2)) + G(pos_ac3(p3,1),pos_ac3(p3,1))*Edio_re(pos_ac3(p3,1)) - B(pos_ac3(p3,1),pos_ac3(p3,2))*Edio_im(pos_ac3(p3,2)) - B(pos_ac3(p3,1),pos_ac3(p3,1))*Edio_im(pos_ac3(p3,1))) + ...
                    Edio_im(pos_ac3(p3,1)).*(B(pos_ac3(p3,1),pos_ac3(p3,2))*Edio_re(pos_ac3(p3,2)) + B(pos_ac3(p3,1),pos_ac3(p3,1))*Edio_re(pos_ac3(p3,1)) + G(pos_ac3(p3,1),pos_ac3(p3,2))*Edio_im(pos_ac3(p3,2)) + G(pos_ac3(p3,1),pos_ac3(p3,1))*Edio_im(pos_ac3(p3,1)))) + ...
                   (-beta*Edio_re(pos_ac3(p3,2)) + beta*Edio_re(pos_ac3(p3,1)) + gamma*Edio_im(pos_ac3(p3,2)) - gamma*Edio_im(pos_ac3(p3,1)));

               
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
    
            alpha =  ( G(pos_dc3(p,1),pos_dc3(p,2))* Edio_re(pos_dc3(p,2)))^2 - 4 * G(pos_dc3(p,1),pos_dc3(p,1)) .* (P_0pn);

            for o = 1:length(alpha)
                if alpha(o) < 1e-12
                    alpha(o) = 1;
                    disp('Warning: alpha is wrong')
                end
            end

%           Ssysm =   Atot*E_star.*conj(Atot*(G+1i*B)*E_star)
            
            A1_re = -1./sqrt(alpha).* ( G(pos_ac3(p3,1) ,pos_ac3(p3,2) )*Edio_re(pos_ac3(p3,1) ) + B(pos_ac3(p3,1) ,pos_ac3(p3,2) )*Edio_im(pos_ac3(p3,1) ) - diag(beta));
            A1_im = -1./sqrt(alpha).* (-B(pos_ac3(p3,1) ,pos_ac3(p3,2) )*Edio_re(pos_ac3(p3,1) ) + G(pos_ac3(p3,1) ,pos_ac3(p3,2) )*Edio_im(pos_ac3(p3,1) ) + diag(gamma));

            A2_re = -1./sqrt(alpha).* ( G(pos_ac3(p3,1) ,pos_ac3(p3,2) )*Edio_re(pos_ac3(p3,2) ) + 2*G(pos_ac3(p3,1) ,pos_ac3(p3,1) )*Edio_re(pos_ac3(p3,1) ) - B(pos_ac3(p3,1) ,pos_ac3(p3,2) )*Edio_im(pos_ac3(p3,2) ) + diag(beta));
            A2_im = -1./sqrt(alpha).* ( B(pos_ac3(p3,1) ,pos_ac3(p3,2) )*Edio_re(pos_ac3(p3,2) ) + G(pos_ac3(p3,1) ,pos_ac3(p3,2) )*Edio_im(pos_ac3(p3,2) ) + 2*G(pos_ac3(p3,1) ,pos_ac3(p3,1) )*Edio_im(pos_ac3(p3,1) ) - diag(gamma));


            J_AFErR(pos_ac3(p3(2),1),pos_ac3(p3,2) ) = real(1/3*inv(A)*complex(A1_re,A1_im))';
            J_AFErR(pos_ac3(p3(2),1),pos_ac3(p3,1) ) = real(1/3*inv(A)*complex(A2_re,A2_im))';
            J_AFErR(pos_ac3(p3(2),1),pos_dc3(p,2) ) = -G(pos_dc3(p,1) ,pos_dc3(p,2) )/(2*G(pos_dc3(p,1) ,pos_dc3(p,1) )) + G(pos_dc3(p,1) ,pos_dc3(p,2) )^2 * Edio_re(pos_dc3(p,2) ) / ( 2* G(pos_dc3(p,1) ,pos_dc3(p,1) ) * sqrt(alpha(2)) );

            J_AFErX(pos_ac3(p3(2),1),pos_ac3(p3,2) ) = imag(1/3*inv(A)*complex(A1_re,A1_im))';
            J_AFErX(pos_ac3(p3(2),1),pos_ac3(p3,1) ) = imag(1/3*inv(A)*complex(A2_re,A2_im))';

            J_AFErR(pos_ac3(p3([1,3]),1),: ) = J_EdiorR(pos_ac3(p3([1,3]),1),: );
            J_AFErX(pos_ac3(p3([1,3]),1),: ) = J_EdiorX(pos_ac3(p3([1,3]),1),: );

        %% Q

        
            J_AFEiR(pos_ac3(p3([2]),1),: ) = J_QpR(pos_ac3(p3([2]),1),: );
            J_AFEiX(pos_ac3(p3([2]),1),: ) = J_QpX(pos_ac3(p3([2]),1),: );

            J_AFEiR(pos_ac3(p3([1,3]),1),: ) = J_EdioiR(pos_ac3(p3([1,3]),1),: );
            J_AFEiX(pos_ac3(p3([1,3]),1),: ) = J_EdioiX(pos_ac3(p3([1,3]),1),: );


        elseif AFE_type(p) == "PQ"

        %% P  

            J_AFErR(pos_ac3(p3([2]),1),: ) = J_PpR(pos_ac3(p3([2]),1),: );
            J_AFErX(pos_ac3(p3([2]),1),: ) = J_PpX(pos_ac3(p3([2]),1),: );
            
            J_AFErR(pos_ac3(p3(2),1),pos_ac3(p3,1) ) = J_AFErR(pos_ac3(p3(2),1),pos_ac3(p3,1) ) + 2*R*E_re(pos_ac3(p3,1))' - 2*R*E_re(pos_ac3(p3,2))' - 2*X*E_im(pos_ac3(p3,1))' + 2*X*E_im(pos_ac3(p3,2))';
            J_AFErR(pos_ac3(p3(2),1),pos_ac3(p3,2) ) = J_AFErR(pos_ac3(p3(2),1),pos_ac3(p3,2) ) - 2*R*E_re(pos_ac3(p3,1))' + 2*R*E_re(pos_ac3(p3,2))' + 2*X*E_im(pos_ac3(p3,1))' - 2*X*E_im(pos_ac3(p3,2))';
            J_AFErX(pos_ac3(p3(2),1),pos_ac3(p3,1) ) = J_AFErX(pos_ac3(p3(2),1),pos_ac3(p3,1) ) - 2*X*E_re(pos_ac3(p3,1))' + 2*X*E_re(pos_ac3(p3,2))' - 2*R*E_im(pos_ac3(p3,1))' + 2*R*E_im(pos_ac3(p3,2))';
            J_AFErX(pos_ac3(p3(2),1),pos_ac3(p3,2) ) = J_AFErX(pos_ac3(p3(2),1),pos_ac3(p3,2) ) + 2*X*E_re(pos_ac3(p3,1))' - 2*X*E_re(pos_ac3(p3,2))' + 2*R*E_im(pos_ac3(p3,1))' - 2*R*E_im(pos_ac3(p3,2))';
           

            J_AFErR(pos_ac3(p3([1,3]),1),: ) = J_EdiorR(pos_ac3(p3([1,3]),1),: );
            J_AFErX(pos_ac3(p3([1,3]),1),: ) = J_EdiorX(pos_ac3(p3([1,3]),1),: );

            
            
        %% Q

            J_AFEiR(pos_ac3(p3([2]),1),: ) = J_QpR(pos_ac3(p3([2]),1),: );
            J_AFEiX(pos_ac3(p3([2]),1),: ) = J_QpX(pos_ac3(p3([2]),1),: );

            J_AFEiR(pos_ac3(p3([1,3]),1),: ) = J_EdioiR(pos_ac3(p3([1,3]),1),: );
            J_AFEiX(pos_ac3(p3([1,3]),1),: ) = J_EdioiX(pos_ac3(p3([1,3]),1),: );
        end
    end

end