function new_state = next_state(state,n_states)
    new_state = mod(state,n_states)+1;
end