function U = EBE_FEM_POISSON_1D(ND,NE,Node,Elem,ua,ub)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         Matlab code for 1D FEM for
%         
%         -u'' = f(x), a<= x <=b, u(a)=ua, u(b)=ub
%
%      Input:  ND - Number of total nodal points
%              NE - Number of total elements
%              Node - Nodal points
%              Elem - Element information
%              ua = Dirichlet BC at x = a
%              ub = Dirichlet BC at x = b
%      Output: U = FEM solution at nodal points
%
%      Function needed: f(x)
%
%      Matlab functions used
%
%      hat1(x,xL,Xr): hat function in [xL,xR] that is 1
%      at xL, and 0 at xR.
%  
%      hat2(x,xL,xR): hat function in [xL,xR] that is 0
%      at xL, and 1 at xR
%
%      int_hat_f(xL,xR): contruction to the load vector from hat1
%
%      int_hat_f(xL,xR): contruction to the load vector from hat 2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


format long;

% step 2: calculate local stiffness matrices/load vectors

LA = zeros(NE,2,2); % initialize a 3D matrix "LA"
LF = zeros(NE,2);

for j = 1:NE
    xL = Node(Elem(j,1));
    xR = Node(Elem(j,2));
    h  = xR - xL;
    
    LA(j,1,1) = 1.0/h;
    LA(j,1,2) = -1.0/h;
    LA(j,2,1) = -1.0/h;
    LA(j,2,2) = 1.0/h;
    
    LF(j,1) = int_hat1_f(xL,xR);
    LF(j,2) = int_hat2_f(xL,xR);
    
end

%Step 3: Assemble global stifffness matrix A and global load vector F


A = zeros(ND,ND);
F = zeros(ND,1);

for j = 1:NE
    nL = Elem(j,1);
    nR = Elem(j,2);
    
    A(nL,nL) = A(nL,nL) + LA(j,1,1);
    A(nL,nR) = A(nL,nR) + LA(j,1,2);
    A(nR,nL) = A(nR,nL) + LA(j,2,1);
    A(nR,nR) = A(nR,nR) + LA(j,2,2);
    
    F(nL) = F(nL) + LF(j,1);
    F(nR) = F(nR) + LF(j,2);
    
end

%Step 4: Impose boundary conditions

for j = 1:ND
    A(1,j) = 0.0;
    A(ND,j) = 0.0;
    
end

A(1,1) = 1.0;
A(ND,ND) = 1.0;

F(1)  = ua;
F(ND) = ub;

%Step 5: Solve the equation by a direct method
U = A\F;

return




