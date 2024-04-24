clear; close all;
%% Full PDE Sim - Electric
%Electrostatic field for heat equation
modelElectric = createpde();
figure
gm = fegeometry("fuse.stl");
pdegplot(gm,"EdgeLabels","on");
importGeometry(modelElectric,"fuse.stl"); % geometryFromEdges for 2-D
Emesh = generateMesh(modelElectric,"Hmax",0.03); % generate mesh 

Conductivity = 2.633E+04; %S/mm

%Specify Equation
specifyCoefficients(modelElectric,'m',0,...
                          'd',0,...
                          'c',Conductivity,...
                          'a',0,...
                          'f',0);

% dirichlet boundary condition for E-domain
% electromagneticProperties(modelElectric,"Conductivity",2.633E+07);
applyBoundaryCondition(modelElectric,"dirichlet","Edge",3,"u",0);

%Neuman boundary condition for E-domain
current = 130 ;%amp
thickness = 1; %mm
boundary_length = 10; %mm
applyBoundaryCondition(modelElectric,"neumann","Edge",6,"q",0,"g",current/(thickness*boundary_length));
% electromagneticBC(modelElectric,"Edge",6, ...current/(thickness*boundary_length)
%     "SurfaceCurrentDensity",current/(thickness*boundary_length)); %surface current density
 
Eresults = solvepde(modelElectric);

pdeplot(modelElectric,"XYData",((Conductivity*Eresults.XGradients).^2+(Conductivity*Eresults.YGradients).^2), ...
              "Contour","on", ...
              "FlowData",[Eresults.XGradients,Eresults.YGradients])
axis equal

%% Full PDE Sim - Heat Set Up
rho_0 = 1/Conductivity; %ohm-mm
alpha = 0.004; % 1/K
f = @(region,state)heatGen_lin(region,state,Eresults,rho_0,thickness, alpha);

% Thermal simulation (Joule heating)
modelThermal = createpde();
% ... (Define geometry, coefficients, and boundary conditions for thermal simulation)
importGeometry(modelThermal,"fuse.stl"); % geometryFromEdges for 2-D
thermal_mesh = generateMesh(modelThermal, "Hmax", 0.5); % generate mesh 

% Specify Joule heating term using electric potential obtained
thermalConductivity = 0.157; % thermal conductivity - W/mmK
rho = 2.7e-6; % mass density - kg/mm^3
SpecificHeat = 978; % J/KgK
enthalpyConstant = rho*SpecificHeat;
specifyCoefficients(modelThermal, "m", 0, "d", enthalpyConstant, "c", thermalConductivity, "a", 0, "f", f);

% Set boundary conditions for thermal simulation
% thermalBC(modelThermal, "Edge", 3, "Temperature", 280); 
applyBoundaryCondition(modelThermal,"dirichlet","Edge",3,"u",60);% Set temperature at the boundary
applyBoundaryCondition(modelThermal,"dirichlet","Edge",6,"u",80);% Set temperature at the boundary

% Generate mesh and solve the thermal PDE
tlist = linspace(0, 1.5, 100); % Define the time vector
setInitialConditions(modelThermal,70);



%% Add source vector data
F = [];
modes_accumulator = [];
for k = 1:100
    state.time=tlist(k);
    state.u = results_thermal.NodalSolution(:,k)';
    MATs = assembleFEMatrices(modelThermal, "nullspace", state);
    F = [F MATs.Fc];
end

uc = MATs.B'*results_thermal.NodalSolution;
norm_uc = norm(uc);
norm_F = norm(F);

modes_accumulator = [modes_accumulator  uc*1000/norm_uc F*1000/norm_F];


%% Conduction only modes
specifyCoefficients(modelThermal, "m", 0, "d", enthalpyConstant, "c", thermalConductivity, "a", 0, "f", 0);
results_thermal = solvepde(modelThermal,tlist);


modes_accumulator = [modes_accumulator MATs.B'*results_thermal.NodalSolution];

% %Plot the temperature distribution
% for i = 1:length(tlist)
%     pdeplot(modelThermal, "XYData", results_thermal.NodalSolution(:,i), "Contour", "on","ColorMap","Hot");
%     title("FOM: Joule Heating 2D Model with Electric Potential Input - PDE Toolbox");
%     axis equal
%     pause(1/10000); % Add a pause if you want to visualize the results in real-time
% end

%% PCA and ROM Modes extraction

X = modes_accumulator;
% X = X - mean(X);
[U, S, ~] = svd(X,"econ");
U_f = FE_Matrices.B*SVD_Matrices.Uc;
for i = 1:10
    subplot(2,5,i)
    pdeplot(modelThermal, "XYData", U_f(:,i))
    % axis equal
    title(strcat("SVD: mode ", string(i)))
end

Ur = U(:, 1:4);


y = zeros(length(S),1);
for i = 1:length(S)
    y(i) = S(i,i);
end

%% Extracting Matrices and Simulating using ODE 45
MATs = assembleFEMatrices(modelThermal, state);

[B,Or] = pdenullorth(MATs.H);
ud = Or*((MATs.H*Or\MATs.R)); % Vector with known value of the constraint DoF.
Kc = B'*(MATs.K + MATs.A + MATs.Q)*B;
Mc = B'*MATs.M*B;
Fc = B'*((MATs.F + MATs.G)-(MATs.K + MATs.A + MATs.Q)*ud);
Shape = size(Fc);

% fun = @(t, u) ((Fc - Kc*u)) ;
% u0 = 70*ones(Shape);
% 
% options = odeset('Mass',Mc ,'AbsTol',1e-3 ,'RelTol',1e-2 ,'Stats','on');
% [t,u] = ode45(fun, tlist, u0, options);
% 
% u=u';
% u = B*u+ud;
% 
% for i = 1:length(tlist)
%     pdeplot(modelThermal, "XYData", u(:,i), "Contour", "on","ColorMap","Hot");
%     title("FOM: Joule Heating 2D Model with Electric Potential Input- ODE45");
%     axis equal
%     pause(1/10000); % Add a pause if you want to visualize the results in real-time
% end

%% ROM Generation
Urc = B'*Ur;
% urd = Urc'*ud; % Vector with known value of the constraint DoF.
Krc = Urc'*Kc*Urc;
Mrc = Urc'*Mc*Urc;
Frc = Urc'*Fc;
ur0 = Urc'*u0;

funr = @(t, u) ((Frc - Krc*u)) ;
options = odeset('Mass',Mrc ,'AbsTol',1e-3 ,'RelTol',1e-3 ,'Stats','on');
[t,ur] = ode45(funr, tlist, ur0, options);

ur=ur';
ur = Urc*ur;
ur = B*ur+ud;

for i = 1:length(tlist)
    pdeplot(modelThermal, "XYData", ur(:,i), "Contour", "on","ColorMap","Hot");
    title("ROM: Joule Heating 2D Model with Electric Potential Input- ODE45");
    axis equal
    pause(1/10000); % Add a pause if you want to visualize the results in real-time
end

%% Error Quantification
err = u - ur;
err_pde = results_thermal.NodalSolution - u;

for i = 1:length(tlist)
    pdeplot(modelThermal, "XYData", err(:,i), "Contour", "on");
    title("ROM to FOM error - Degrees Celcius");
    axis equal
    pause(1/10000); % Add a pause if you want to visualize the results in real-time
end

for i = 1:length(tlist)
    pdeplot(modelThermal, "XYData", err_pde(:,i), "Contour", "on");
    title("PDE to ODE45 error - Degrees Celcius");
    axis equal
    pause(1/100); % Add a pause if you want to visualize the results in real-time
end

%% 
sensors = findNodes(modelThermal.Mesh,"nearest", [0,0,0; 0, 5, 10]);
projection_mat = zeros(length(sensors), size(results_thermal.NodalSolution,1));
for i = 1:length(sensors)
    projection_mat(sensors(i)) = 1;
end
projection_mat = sparse(projection_mat);

temp_sensors = projection_mat*results_thermal.NodalSolution;
