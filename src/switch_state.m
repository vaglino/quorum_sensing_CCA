function states = switch_state(states,n_states,positive)
    % switch colonies to their next state next state is determined with
    % cyclical numbers. i.e. for 3 states, the transition can only be 
    % 1->2->3->1... 
    % this is mathematically equivalent to (i-1 modulo 3) + 1
    
    new_states = mod(states(positive),n_states)+1;
%     new_states = next_state(states(positive),n_states);
    states(positive) = new_states;
end