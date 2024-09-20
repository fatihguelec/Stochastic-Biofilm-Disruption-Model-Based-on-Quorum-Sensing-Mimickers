function [j, tau] = gillespie_direct(a, N_reac)
        % Gillespie Direct Algorithm
        s_1 = rand; s_2 = rand; % Generate random numbers from U(0,1)
        a0 = sum(a,'all'); 
        tau = (1/a0)*log(1/s_1); %calculate when the reaction takes place
        
        % Determine which reaction to take place
        for j = 1:N_reac
            if sum(a(1:j)) >= s_2*a0 && sum(a(1:j-1)) < s_2*a0
                break;
            end
        end
end