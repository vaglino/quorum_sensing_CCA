clear all
% close all
% % Determine where your m-file's folder is.
% folder = fileparts(which(mfilename)); 
% % Add that folder plus all subfolders to the path.
% addpath(genpath(folder));

%% parameters
g = struct;  % struture with geometrical parameters
g.nx = 100;                           % Number of steps in space(x)
g.ny = 100;                           % Number of steps in space(y)
g.dx=2/(g.nx-1);                     % Width of space step(x)
g.dy=2/(g.ny-1);                     % Width of space step(y)
g.x=0:g.dx:2;                        % Range of x(0,2) and specifying the grid points
g.y=0:g.dy:2;                        % Range of y(0,2) and specifying the grid points

p = struct;  % other parameters
p.nt = 100;                          % Number of time steps 
p.dt = 0.01;                         % Width of each time step
p.n_states = 4;                      % number of states
p.diff = 0.1;                         % Diffusion coefficient/viscocity
p.gamma = 0.01;                      % decay rate
p.production_rate = 0.2;             % rates at which bacteria produces AI
p.thresh = 0.3;                      % AI sensing threshold

u=zeros(g.nx,g.ny,p.n_states);       % Prealocating u (AI concentration field)
% un=zeros(g.nx,g.ny);                 % Preallocating un

Dir.UW=0;                            % x=0 Dirichlet B.C 
Dir.UE=0;                            % x=L Dirichlet B.C 
Dir.US=0;                            % y=0 Dirichlet B.C 
Dir.UN=0;                            % y=L Dirichlet B.C 
Neu.UW=0;                            % x=0 Neumann B.C (du/dn=UnW)
Neu.UE=0;                            % x=L Neumann B.C (du/dn=UnE)
Neu.US=0;                            % y=0 Neumann B.C (du/dn=UnS)
Neu.UN=0;                            % y=L Neumann B.C (du/dn=UnN)

BC_type = "Neumann";
% BC = "Dirichlet"

%% initial conditions
% randomly initiated colonies on grid
% Initial Conditions
% for i=1:g.nx
%     for j=1:g.ny
%         if ((1<=g.y(j))&&(g.y(j)<=1.5)&&(1<=g.x(i))&&(g.x(i)<=1.5))
%             u(i,j)=2;
%         else
%             u(i,j)=0;
%         end
%     end
% end
% u(10,10) = 1;

% u=zeros(g.nx,g.ny);                  %Preallocating u
                  %Preallocating u

% sour = double(rand(size(u))>=0.9)*0.05;
% sour2 = zeros(size(u));

% sour = zeros(size(u));
% sour(10,10) = 0.1;
% sour(3,3) = 0.1;
% sour(20:21,30:31) = 2;
% sour(15:17,70:72) = 2;

% colonies = sour ~= 0;


colonies = rand(g.nx,g.ny) >= 0.9;
% states = randi([1, p.n_states],g.nx,g.ny) .* colonies;
states = randi([1, 1],g.nx,g.ny) .* colonies;

% define boundary conditions
bc = define_bc(g,p,BC_type,Neu);
D = coeff_matrix(g,p,BC_type);

run_simulation(u,D,bc,states,g,p,BC_type,Neu);

% maxu = [];
% for it=0:p.nt
%     u(sour~=0) = u(sour~=0) + sour(sour~=0);
%     plot_field(g.x,g.y,u,p,it)
%     u = diffuse(u,D,bc,g,p,BC_type,Neu);
%     positive_colonies = colonies_above_threshold(colonies,u,thresh);
%     sour = end_sources(sour,positive_colonies);
% %     plot_field(g.x,g.y,u,p,it)
%     maxu = [maxu; max(max(u))];
% end

%% simulate system

function run_simulation(u,D,bc,states,g,p,BC_type,bounds)

for t=0:p.nt % at every time step
    
    for i=1:size(u,3) % for every QS autoinducer
        ui = u(:,:,i);
        % update QSAI concentration field with colonies producing QSAI
        ui(states==i) = ui(states==i) + p.production_rate;
        
        % plot concentration field for QSAI
        subplot(3,2,i)
        plot_field(g.x,g.y,ui,p,t,i)
        
        % take one diffusion step
        ui = diffuse(ui,D,bc,g,p,BC_type,bounds);
   
        u(:,:,i) = ui;
    end
    % figure out if any colony senses QSAI concentration above
    % threshold
    positive_colonies = colonies_above_threshold(states,u,p.thresh);
    % update colonies to new states
    states = switch_state(states,p.n_states,positive_colonies);
    
    subplot(3,2,i+1)
    plot_superimposed(g.x,g.y,states,g,p,t,i)
end
end

%% B.C vector
function [bc] = define_bc(g,p,BC_type,edge)
    bc=zeros(g.nx-2,g.ny-2);
    if BC_type == "Dirichlet"
        bc(1,:)=edge.UW/g.dx^2; bc(g.nx-2,:)=edge.UE/g.dx^2;  %Dirichlet B.Cs
        bc(:,1)=edge.US/g.dy^2; bc(:,g.ny-2)=edge.UN/g.dy^2;  %Dirichlet B.Cs
    elseif BC_type == "Neumann"
        bc(1,:)=-edge.UW/g.dx; bc(g.nx-2,:)=edge.UE/g.dx;  %Neumann B.Cs
        bc(:,1)=-edge.US/g.dy; bc(:,g.nx-2)=edge.UN/g.dy;  %Neumann B.Cs
    end
    %B.Cs at the corners:
    bc(1,1)=edge.UW/g.dx^2+edge.US/g.dy^2; 
    bc(g.nx-2,1)=edge.UE/g.dx^2+edge.US/g.dy^2;
    bc(1,g.ny-2)=edge.UW/g.dx^2+edge.UN/g.dy^2; 
    bc(g.nx-2,g.ny-2)=edge.UE/g.dx^2+edge.UN/g.dy^2;
    
    bc=p.diff*p.dt*bc;
end

%% Calculating the coefficient matrix for the implicit scheme
function D = coeff_matrix(g,p,BC_type)
    Ex=sparse(2:g.nx-2,1:g.nx-3,1,g.nx-2,g.nx-2);
    Ax=Ex+Ex'-2*speye(g.nx-2);        %Dirichlet B.Cs
    if BC_type == "Neumann"  
        Ax(1,1)=-1; Ax(g.nx-2,g.nx-2)=-1;  %Neumann B.Cs
    end
    Ey=sparse(2:g.ny-2,1:g.ny-3,1,g.ny-2,g.ny-2);
    Ay=Ey+Ey'-2*speye(g.ny-2);        %Dirichlet B.Cs
    if BC_type == "Neumann"  
        Ay(1,1)=-1; Ay(g.ny-2,g.ny-2)=-1;  %Neumann B.Cs
    end
    A=kron(Ay/g.dy^2,speye(g.nx-2))+kron(speye(g.ny-2),Ax/g.dx^2);
    D=speye((g.nx-2)*(g.ny-2))-p.diff*p.dt*A;
end

function u = diffuse(u,D,bc,g,p,BC_type,edge)
    
    un=u;
    U=un;U(1,:)=[];U(end,:)=[];U(:,1)=[];U(:,end)=[];
    U=reshape(U+bc,[],1);
    
    % diffusion step
    U = D\U;
    U = reshape(U,g.nx-2,g.ny-2);
    u(2:g.nx-1,2:g.ny-1) = U;
    
    %Boundary conditions
    if BC_type == "Dirichlet"
    %Dirichlet:
        u(1,:)   = edge.UW;
        u(end,:) = edge.UE;
        u(:,1)   = edge.US;
        u(:,end) = edge.UN;
    %Neumann:
    elseif BC_type == "Neumann"
        u(1,:)   = u(2,:)     - edge.UW*g.dx;
        u(end,:) = u(end-1,:) + edge.UE*g.dx;
        u(:,1)   = u(:,2)     - edge.US*g.dy;
        u(:,end) = u(:,end-1) + edge.UN*g.dy;
    end
    u = u - p.gamma*u;
end

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
    
function sour = end_sources(sour,positive_colonies)
    sour(positive_colonies) = 0;
end

function states = switch_state(states,n_states,positive)
    % switch colonies to their next state next state is determined with
    % cyclical numbers. i.e. for 3 states, the transition can only be 
    % 1->2->3->1... 
    % this is mathematically equivalent to (i-1 modulo 3) + 1
    
    new_states = mod(states(positive),n_states)+1; 
    states(positive) = new_states;
end