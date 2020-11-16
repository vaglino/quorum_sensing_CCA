function positive_colonies = colonies_above_threshold(states,u,p)

    positive_colonies = zeros(size(states));
    n_states = size(u,3);
    for i=1:n_states
        % find which colonies are at state i
        i_colonies = states == i; 
        % calculate AI_i concentration each colony at state i senses
        next_i = next_state(i,p.n_states);
        
        % check concentration of next AI. if colony is i, check concentration
        % of (i+1)
        u_at_colonies = u(:,:,next_i) .* i_colonies;
        % check concentration of same AI. if colony is i, check concentration
        % of i
%         u_at_colonies = u(:,:,i) .* i_colonies;  

        positive_colonies_i =  u_at_colonies >= p.thresh;
        positive_colonies = positive_colonies + positive_colonies_i;
    end
    positive_colonies = logical(positive_colonies);
end