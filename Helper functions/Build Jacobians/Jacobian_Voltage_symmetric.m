function [J_EdiorR, J_EdiorX, J_EdioiR, J_EdioiX] = Jacobian_Voltage_symmetric(Grid_para, Filter_para, J_EdiorR, J_EdiorX, J_EdioiR, J_EdioiX)

    n_dc = Grid_para.n_dc;
    n_ac = Grid_para.n_ac;
    n_nodes = Grid_para.n_nodes;
    n_ph = Grid_para.n_ph;
    R = Filter_para.R;
    X = Filter_para.X;
    
    pos_ac3 = Grid_para.pos_ac3;
    pos_dc3 = Grid_para.pos_dc3;
    
    alp = exp(2*pi/3*1i);
    A = 1/3*[1 1     1; 
             1 alp   alp^2; 
             1 alp^2 alp];
         
    gamma = ((1-R) - 1i*X) * A;
    delta = (R + 1i*X) * A;
    G1 = real(gamma);
    G2 = -imag(gamma);
    G3 = imag(gamma);
    G4 = real(gamma);
    
    D1 = real(delta);
    D2 = -imag(delta);
    D3 = imag(delta);
    D4 = real(delta);
    
    A1 = 1/3*[1 1 1; 1 real(alp)   -imag(alp);   1 real(alp^2) -imag(alp^2)];
    A2 = 1/3*[0 0 0; 0 real(alp^2) -imag(alp^2); 0 real(alp)   -imag(alp)];
    A3 = 1/3*[0 0 0; 1 imag(alp)    real(alp);   1 imag(alp^2)  real(alp^2)];
    A4 = 1/3*[1 1 1; 0 imag(alp^2)  real(alp^2); 0 imag(alp)    real(alp)];   

    A1 = real(A);
    A2 = -imag(A);
    A3 = imag(A);
    A4 = real(A);
    
    for i=1:n_ac
            i3 = polyphase_indices(i,n_ph);

            
            J_EdiorR(i3,i3) = A1;
            J_EdiorX(i3,i3) = A2;
            J_EdioiR(i3,i3) = A3;
            J_EdioiX(i3,i3) = A4;
    end

%     for i = 1:size(pos_ac3,1)/Grid_para.n_ph
%             i3 = polyphase_indices(i,n_ph);
%         
%             J_EdiorR(pos_ac3(i3,1), pos_ac3(i3,1)) = G1;
%             J_EdiorX(pos_ac3(i3,1), pos_ac3(i3,1)) = G2;
%             J_EdioiR(pos_ac3(i3,1), pos_ac3(i3,1)) = G3;
%             J_EdioiX(pos_ac3(i3,1), pos_ac3(i3,1)) = G4;
%             
%             J_EdiorR(pos_ac3(i3,1), pos_ac3(i3,2)) = D1;
%             J_EdiorX(pos_ac3(i3,1), pos_ac3(i3,2)) = D2;
%             J_EdioiR(pos_ac3(i3,1), pos_ac3(i3,2)) = D3;
%             J_EdioiX(pos_ac3(i3,1), pos_ac3(i3,2)) = D4;
%     end
    
end

