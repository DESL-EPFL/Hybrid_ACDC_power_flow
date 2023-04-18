function [J_PpR, J_PpX, J_QpR, J_QpX] = Jacobian_Powers_symmetric(E_re,E_im, Grid_para, J_PpR, J_PpX, J_QpR, J_QpX)

    n_dc = Grid_para.n_dc;
    n_ac = Grid_para.n_ac;
    n_nodes = Grid_para.n_nodes;
    n_ph = Grid_para.n_ph;
    G = Grid_para.G;
    B = Grid_para.B;
        
    alp = exp(2*pi/3*1i);
    A = 1/3*[1 1     1; 
             1 alp   alp^2; 
             1 alp^2 alp];

    ACell =  repmat({A}, 1, n_ac);
    ICell =  repmat({1}, 1, n_dc); 
    Atot = blkdiag(ACell{:},ICell{:});

    Edio_re = real(Atot*(complex(E_re,E_im)));
    Edio_im = imag(Atot*(complex(E_re,E_im)));

    for i=1:n_nodes

        % Diagonal elements (terms outside the sum)
        J_PpR(i,i) = 2*G(i,i)*Edio_re(i);
        J_PpX(i,i) = 2*G(i,i)*Edio_im(i);
        J_QpR(i,i) = -2*B(i,i)*Edio_re(i);
        J_QpX(i,i) = -2*B(i,i)*Edio_im(i);

        for j=1:n_nodes
            if(j~=i)

                % Diagonal elements (terms inside the sum)
                J_PpR(i,i) = J_PpR(i,i) + G(i,j)*Edio_re(j) - B(i,j)*Edio_im(j);
                J_PpX(i,i) = J_PpX(i,i) + B(i,j)*Edio_re(j) + G(i,j)*Edio_im(j);
                J_QpR(i,i) = J_QpR(i,i) - B(i,j)*Edio_re(j) - G(i,j)*Edio_im(j);
                J_QpX(i,i) = J_QpX(i,i) + G(i,j)*Edio_re(j) - B(i,j)*Edio_im(j);

                % Offdiagonal elements
                J_PpR(i,j) =  G(i,j)*Edio_re(i)+B(i,j)*Edio_im(i);
                J_PpX(i,j) = -B(i,j)*Edio_re(i)+G(i,j)*Edio_im(i);
                J_QpR(i,j) = -B(i,j)*Edio_re(i)+G(i,j)*Edio_im(i);
                J_QpX(i,j) = -G(i,j)*Edio_re(i)-B(i,j)*Edio_im(i);
            end 
        end
    end

    J_SpR = complex(J_PpR,J_QpR);
    J_SpX = complex(J_PpX,J_QpX);

    J_SpR2 = zeros(size(J_PpR));
    J_SpX2 = zeros(size(J_PpR));
    for i = 1:n_ac
        for j = 1:n_ac
            i3 = polyphase_indices(i,n_ph);
            j3 = polyphase_indices(j,n_ph);

            J_SpR2(i3,j3) = 1/3*inv(A).*diag(J_SpR(i3,j3));
            J_SpX2(i3,j3) = 1/3*inv(A).*diag(J_SpX(i3,j3));
        end
    end

    J_PpR = real(J_SpR2);
    J_QpR = imag(J_SpR2);

    J_PpX = real(J_SpX2);
    J_QpX = imag(J_SpX2);

end
