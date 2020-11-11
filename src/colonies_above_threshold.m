function positive_colonies = colonies_above_threshold(states,u,thresh)
    positive_colonies = zeros(size(states));
    n_states = size(u,3);
    for i=1:n_states
        % find which colonies are at state i
        i_colonies = states == i; 
        % calculate AI_i concentration each colony at state i senses
        u_at_colonies = u(:,:,i) .* i_colonies;  
        positive_colonies_i =  u_at_colonies >= thresh;
        positive_colonies = positive_colonies + positive_colonies_i;
    end
    positive_colonies = logical(positive_colonies);
end