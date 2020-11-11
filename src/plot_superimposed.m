function plot_superimposed(x,y,states,g,p,t,i)
    colors = {'b','r','g','m'};
  
    for i=1:p.n_states
        [row,col,v] = find(states==i);
        h = scatter(row*g.dy,col*g.dx,10,colors{i},'filled','square');
        hold on
    end%
    axis ([0 2 0 2 0 2])
    view(2)
    pbaspect([1 1 1])
    title({['states matrix'];['time (\itt) = ',num2str(t*p.dt)]})
    xlabel('X \rightarrow')
    ylabel('{\leftarrow} Y')
    drawnow; 
    refreshdata(h)
    hold off
end