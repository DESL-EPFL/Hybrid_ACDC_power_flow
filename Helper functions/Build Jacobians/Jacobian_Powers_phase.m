 function [J_PR, J_PX, J_QR, J_QX, J_ER, J_EX] = Jacobian_Powers_phase(E_re,E_im,Grid_para, J_PR, J_PX, J_QR, J_QX, J_ER, J_EX)

        n_dc = Grid_para.n_dc;
        n_ac = Grid_para.n_ac;
        n_nodes = Grid_para.n_nodes;
        n_ph = Grid_para.n_ph;
        G = Grid_para.G;
        B = Grid_para.B;
        
        for i=1:n_nodes

            % Diagonal elements (terms outside the sum)
            J_PR(i,i) = 2*G(i,i)*E_re(i);
            J_PX(i,i) = 2*G(i,i)*E_im(i);
            J_QR(i,i) = -2*B(i,i)*E_re(i);
            J_QX(i,i) = -2*B(i,i)*E_im(i);
            J_ER(i,i) = 2*E_re(i);
            J_EX(i,i) = 2*E_im(i);

            if i > n_ac*n_ph
                J_ER(i,i) = 1;
                J_EX(i,i) = 0;
            end

            for j=1:n_nodes
                if(j~=i)

                    % Diagonal elements (terms inside the sum)
                    J_PR(i,i) = J_PR(i,i) + G(i,j)*E_re(j) - B(i,j)*E_im(j);
                    J_PX(i,i) = J_PX(i,i) + B(i,j)*E_re(j) + G(i,j)*E_im(j);
                    J_QR(i,i) = J_QR(i,i) - B(i,j)*E_re(j) - G(i,j)*E_im(j);
                    J_QX(i,i) = J_QX(i,i) + G(i,j)*E_re(j) - B(i,j)*E_im(j);

                    % Offdiagonal elements
                    J_PR(i,j) =  G(i,j)*E_re(i)+B(i,j)*E_im(i);
                    J_PX(i,j) = -B(i,j)*E_re(i)+G(i,j)*E_im(i);
                    J_QR(i,j) = -B(i,j)*E_re(i)+G(i,j)*E_im(i);
                    J_QX(i,j) = -G(i,j)*E_re(i)-B(i,j)*E_im(i);
                end 
            end
        end

    end