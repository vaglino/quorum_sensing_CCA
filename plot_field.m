function plot_field(x,y,u,g,p,t,i)
    maps = {'Blues','Reds','Greens','Greys','Purples'};
    map = brewermap(100,maps{i}); 
    map(1,:) = [1,1,1];
    h=surf(x,y,u','EdgeColor','none');       %plotting the field variable
%     shading interp
    colormap(gca,map)
%     axis ([0 2 0 2 0 3])
    axis ([0 g.w 0 g.h 0 3])
    view(2)
    pbaspect([1 1 1])
%     title({['2-D Diffusion with D = ',num2str(p.diff)];['time (\itt) = ',num2str(t*p.dt)]})
    title({['2-D Diffusion of AI i = ',num2str(i)];['time (\itt) = ',num2str(t*p.dt)]})
    xlabel('X \rightarrow')
    ylabel('Y \rightarrow')
    zlabel('Transport property profile (u) \rightarrow')
%     drawnow; 
    
    refreshdata(h)
end