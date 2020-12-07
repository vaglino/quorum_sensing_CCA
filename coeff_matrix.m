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