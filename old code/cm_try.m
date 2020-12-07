% pcolor(3.5*abs(peaks(1024)))
% shading interp
% % Create colormap: 
map = brewermap(8,'GnBu');
% map(1,:) = [1 1 1]; % optionally force first color to white 

[X,Y] = meshgrid(-5:.5:5);
Z = Y.*sin(X) - X.*cos(Y);
s = surf(X,Y,Z,'cm',map)%,'FaceAlpha',0.5)