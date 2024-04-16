function [Data, tlist, FE_Matrices, SVD_Matrices, projection_mat] = genData(curr, BC, sensor_positions)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Data = cell(0); % initialize container for data

%% PDE Sim Set up - Electric
%Electrostatic field for heat equation
modelElectric = createpde();
importGeometry(modelElectric,"fuse.stl"); % geometryFromEdges for 2-D
generateMesh(modelElectric,"Hmax",0.1); % generate mesh 

Conductivity = 2.633E+04; %S/mm

%Specify Equation
specifyCoefficients(modelElectric,'m',0,...
                          'd',0,...
                          'c',Conductivity,...
                          'a',0,...
                          'f',0);

% dirichlet boundary condition for E-domain
applyBoundaryCondition(modelElectric,"dirichlet","Edge",3,"u",0);

%Neuman boundary condition for E-domain

thickness = 1; %mm
boundary_length = 10; %mm

%% PDE Sim - Heat Set Up
rho_0 = 1/Conductivity; %ohm-mm
alpha = 0.004; % 1/K
% Thermal simulation (Joule heating)
modelThermal = createpde();
% ... (Define geometry, coefficients, and boundary conditions for thermal simulation)
importGeometry(modelThermal,"fuse.stl"); % geometryFromEdges for 2-D
generateMesh(modelThermal,"Hmax",0.5); % generate mesh 

% Specify Joule heating term using electric potential obtained
thermalConductivity = 0.157; % thermal conductivity - W/mmK
rho = 2.7e-6; % mass density - kg/mm^3
SpecificHeat = 978; % J/KgK
enthalpyConstant = rho*SpecificHeat;

% Set initial conditions for Thermal PDE
tlist = linspace(0, 2, 100); % Define the time vector
setInitialConditions(modelThermal, 0); % Initializing the Initial conditions to zero since they can be added after the fact
        
modes_accumulator = [];

for i = 1:length(curr)
    for j = 1:size(BC,1)
        %% Electric BCs
        current = curr(i) ;%amp - This is an input from the function
        applyBoundaryCondition(modelElectric,"neumann","Edge",6,"q",0,"g",current/(thickness*boundary_length));
        
        Eresults = solvepde(modelElectric);
        
        %% Thermal BCs
        f = @(region,state)heatGen_lin(region,state,Eresults,rho_0,thickness, alpha);
        specifyCoefficients(modelThermal, "m", 0, "d", enthalpyConstant, "c", thermalConductivity, "a", 0, "f", f);
        
        % Set boundary conditions for thermal simulation
        % thermalBC(modelThermal, "Edge", 3, "Temperature", 280); 
        applyBoundaryCondition(modelThermal,"dirichlet","Edge",3,"u",BC(j,1));% Set temperature at the boundary
        applyBoundaryCondition(modelThermal,"dirichlet","Edge",6,"u",BC(j,2));% Set temperature at the boundary
        
        % Generate mesh and solve the thermal PDE
        results_thermal = solvepde(modelThermal,tlist);
        
        %% Add source vector data
        F = [];
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

        Data = [Data; {[curr(i) BC(j,1) BC(j,2)], uc, F}];
    end 
end


FE_Matrices.Kc = MATs.Kc;
FE_Matrices.Mc = MATs.M;

sensors = findNodes(modelThermal.Mesh,"nearest", sensor_positions);
projection_mat = zeros(length(sensors), size(results_thermal.NodalSolution,1));

for i = 1:length(sensors)
    projection_mat(i,sensors(i)) = 1;
end
projection_mat = sparse(projection_mat);

FE_Matrices.L = projection_mat*MATs.B;

X = modes_accumulator ;
[U, S, ~] = svd(X,"econ");
SVD_Matrices.Uc = U;
SVD_Matrices.Sc = S;



