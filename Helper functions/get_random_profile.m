function S_star = get_random_profile(random_data,S)

    P_min = real(S(1));
    P_max = real(S(2));

    Q_min = imag(S(1));
    Q_max = imag(S(2));
    
    P_star = random_data(1,:) * (-P_min + P_max) +P_min;
    Q_star = random_data(1,:) * (-Q_min + Q_max) +Q_min;
    
    S_star = complex(P_star,Q_star);
    
end
