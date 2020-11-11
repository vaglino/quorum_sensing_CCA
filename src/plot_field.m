function plot_field(x,y,u,p,t,i)
    maps = {'Blues','Reds','Greens','Purples'};
    map = brewermap(100,maps{i}); 
    map(1,:) = [1,1,1];
    h=surf(x,y,u','EdgeColor','none');       %plotting the field variable
%     shading interp
    colormap(gca,map)
    axis ([0 2 0 2 0 2])
    view(2)
    pbaspect([1 1 1])
    title({['2-D Diffusion with D = ',num2str(p.diff)];['time (\itt) = ',num2str(t*p.dt)]})
    xlabel('Spatial co-ordinate (x) \rightarrow')
    ylabel('{\leftarrow} Spatial co-ordinate (y)')
    zlabel('Transport property profile (u) \rightarrow')
%     drawnow; 
    
    refreshdata(h)
end