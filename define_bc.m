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