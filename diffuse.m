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
    u = u - p.gamma*u; % decay with rate gamma
  
end