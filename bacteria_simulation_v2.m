clear all
clf
% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

%% parameters
g = struct;  % struture with geometrical parameters
g.h = 2;
g.w = 2;
g.nx = 300;                           % Number of steps in space(x)
g.ny = 300;                           % Number of steps in space(y)
g.dx=g.w/(g.nx-1);                     % Width of space step(x)
g.dy=g.h/(g.ny-1);                     % Width of space step(y)
g.x=0:g.dx:g.w;                        % Range of x(0,2) and specifying the grid points
g.y=0:g.dy:g.h;                        % Range of y(0,2) and specifying the grid points

p = struct;  % structure with other parameters
p.nt = 100;                          % Number of time steps 
p.dt = 0.02;                         % Width of each time step
p.n_states = 4;                      % number of states
p.diff = 0.05;                         % Diffusion coefficient/viscocity
p.gamma = 0.001;                      % decay rate
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

F = p.diff * p.dt / (g.dx*g.dy)

%% initial conditions

colonies = rand(g.nx,g.ny) >= 0.8; % 10% of squares are colonies 
% states = randi([1, 1],g.nx,g.ny) .* colonies; % start with all colonies in state 1
% uncomment to start with colonies in random state...
states = randi([1, p.n_states],g.nx,g.ny) .* colonies; 

% define boundary conditions
bc = define_bc(g,p,BC_type,Neu); % use neumann boundary conditions
% define coefficient matrix for implicit diffusion solver scheme
D = coeff_matrix(g,p,BC_type);

%% run simulation
run_simulation(u,D,bc,states,g,p,BC_type,Neu);
