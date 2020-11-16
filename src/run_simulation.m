function run_simulation(u,D,bc,states,g,p,BC_type,bounds)

for t=0:p.nt % at every time step
    
    
    subplot(3,2,p.n_states+1)
    plot_superimposed(g.x,g.y,states,g,p,t)
    
    for i=1:p.n_states % for every QS autoinducer
        ui = u(:,:,i);
        % update QSAI concentration field with colonies producing QSAI
        ui(states==i) = ui(states==i) + p.production_rate;
        
        % plot concentration field for QSAI
        subplot(3,2,i)
        plot_field(g.x,g.y,ui,g,p,t,i)
        
        % take one diffusion step
        ui = diffuse(ui,D,bc,g,p,BC_type,bounds);
   
        u(:,:,i) = ui;
    end
    % figure out if any colony senses QSAI concentration above
    % threshold
    positive_colonies = colonies_above_threshold(states,u,p);
    % update colonies to new states
    states = switch_state(states,p.n_states,positive_colonies);
    
    
end
end