clear all       % clear variables and fuctions from memory
close all       % close figures

format long;    % double-precision

% Step 1: Generate the mesh (uniform mesh)

ND = 100;     % number of total nodes
NE = ND-1;   % number of total elements

% Save all node coordinates in the vector "Node"
Node = zeros(ND,1); % Initialize "Node"
h = 10.0/NE;
for j = 1:ND
    Node(j) = (j-1)*h;
end

% Save all elements information in the 2D matrix "Elem"
Elem = zeros(NE,2);  % Initialize "Elem"

for j = 1:NE
    Elem(j,1) = j;
    Elem(j,2) = j+1;
end

%Specify Dirichlet boundary condition:

ua = 40.0;
ub = 200.0;

%Step 2-5: FEM analysis (call EBE_FEM_POISSON_1D)

uAppr = zeros(ND,1);
uAppr = EBE_FEM_POISSON_1D(ND,NE,Node,Elem,ua,ub);

% Step 6: Post-processing, say, plt the numerical solution

figure(1)
plot(Node, uAppr, 'r-o', 'LineWidth', 2);
h_legend = legend('Numerical solution', 2);
set(h_legend,'FontName','Times','FontSize',18,'FontAngle','normal');
xlabel('{\sl{x}}','FontName','Times','FontSize', 18);
ylabel('{\sl{u}}','FontName','Times','FontSize', 18);
hold off
%
%-----------------------------------------------------------
%
%Step 7: Error analysis (if exact solution is available)
%

x = 0:0.01:10.0;  % create vector: x = (0, 0.01, 0.02, 0.03, ..., 1)
K = length(x);   % find length of vector "x"
for i = 1:K
    uExact(i) = -5.0*(x(i).^2) + 66.0.*x(i) + 40.0;
end

for j = 1:ND
    utmp = -5.0*(Node(j).^2) + 66.0.*Node(j) + 40.0;
    error(j) = uAppr(j)-utmp;
end


figure(3)
plot(Node, uAppr,'r-.o','LineWidth', 2);
hold on
plot(x, uExact,'b-','LineWidth', 2);
h_legend = legend('Numerical solution','Exact solution', 2);
set(h_legend,'FontName','Times','FontSize', 18, 'FontAngle','normal');
xlabel('{\sl{x}}','FontName','Times','FontSize', 18);
ylabel('{\sl{u}}','FontName','Times','FontSize', 18);
hold off

figure(5)
plot(Node, error, 'r-o', 'LineWidth', 2);
xlabel('{\sl{x}}','FontName','Times','FontSize', 18);
ylabel('{\sl{Error}}','FontName','Times','FontSize', 18);
hold off


