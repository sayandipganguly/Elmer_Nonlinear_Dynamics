%A Unified code of Elmer and MATLAB**** Written by SayandipGanguly******* This particular example is solved for caltilever beam of square cross-section. Users can use other nonlinear model for the system identification%
clc;clear all;clc;
%define paths
%===================================================================================================================================================================================
 INPUT_ELMERPATH="PROVIDE THE PATH OF ELMER:"; %example of input-'C:\Program Files\Elmer 9.0-Release\bin'
 elmer_path = input(INPUT_ELMERPATH);
 INPUT_PROJECTPATH="PROVIDE THE PATH WHERE GEOMETRIC MODEL IS SAVED:"; %example of input-'D:\Elmer_Axial bar\New folder'
 project_path = input(INPUT_PROJECTPATH);

%create string for line seperator===========================================================================================================================================================================
    lb = sprintf('\n');

% fileid = fopen('geometry.grd','w');

%generate *.mesh files%=============================================================================================================================================================================

%meshfilename='geometry.grd';
%meshpath = ['D:\Elmer_Axial bar\New folder\geometry.grd'];
simulation_command1 = ['ElmerGrid 14 2 geometry -autoclean']; % Replace 14 for Gmsh mesh format by 1,4,5,8,10 respectively for '.grd','.ansys','.inp','.unv','.dat' input file format 
system(simulation_command1);

% To copy mesh files from generated geometry folder==================================================================================================================================   

p1 = [project_path filesep 'geometry'];
% Copying file from respective path
pdest = project_path;
filename1 = 'mesh.boundary';
filename2 = 'mesh.elements';
filename3 = 'mesh.header';
filename4 = 'mesh.nodes';
source1 = fullfile(p1,filename1);
destination1 = fullfile(pdest,filename1);
source2 = fullfile(p1,filename2);
destination2 = fullfile(pdest,filename2);
source3 = fullfile(p1,filename3);
destination3 = fullfile(pdest,filename3);
source4 = fullfile(p1,filename4);
destination4 = fullfile(pdest,filename4);
copyfile(source1,destination1)
copyfile(source2,destination2)
copyfile(source3,destination3)
copyfile(source4,destination4)

%generate *.sif file
%===================================================================================================================================================================================
%define path

    sif_path =  [project_path filesep 'case.sif'];
    fileid = fopen(sif_path,'w');

% write file in .sif
    
    fwrite(fileid,['Header' lb]);
    fwrite(fileid,[' CHECK KEYWORDS Warn' lb]);
    fwrite(fileid,[' Mesh DB "."' lb]);
    fwrite(fileid,[' Include Path "."' lb]);
    fwrite(fileid,[' Results Directory "."' lb]);
    fwrite(fileid,['End' lb]);
    fwrite(fileid, lb);

    fwrite(fileid,['Simulation' lb]);
    fwrite(fileid,[' Max Output Level = 1' lb]);
    fwrite(fileid,[' Coordinate System = Cartesian' lb]);
    fwrite(fileid,[' Coordinate Mapping(3) = 1 2 3' lb]);
    fwrite(fileid,[' Simulation Type = Transient' lb]);
    fwrite(fileid,[' Steady State Max Iterations = 1' lb]);
    fwrite(fileid,[' Output Intervals = 1' lb]);
    fwrite(fileid,[' Timestepping Method = BDF' lb]);
    fwrite(fileid,[' BDF Order = 1' lb]);
    fwrite(fileid,[' Timestep intervals =2000' lb]); %%%Change it as per intended time of analysis
    fwrite(fileid,['  Timestep Sizes = 0.005' lb]);  
    fwrite(fileid,[' Solver Input File = case.sif' lb]);
    fwrite(fileid,[' Post File = case.vtu' lb]); %output file format
    fwrite(fileid,['End' lb]);
    fwrite(fileid, lb);
% 
    fwrite(fileid,['Constants' lb]);
    fwrite(fileid,[' Gravity(4) = 0 -1 0 9.82' lb]);
    fwrite(fileid,[' Stefan Boltzmann = 5.67e-08' lb]);
    fwrite(fileid,[' Permittivity of Vacuum = 8.8542e-12' lb]);
    fwrite(fileid,[' Boltzmann Constant = 1.3807e-23' lb]);
    fwrite(fileid,[' Unit Charge = 1.602e-19' lb]);
    fwrite(fileid,['End' lb]);
    fwrite(fileid, lb);
%Define bodies or part as per the geometric modeling 
    fwrite(fileid,['Body 1' lb]);
    fwrite(fileid,[' Target Bodies(1) = 1' lb]);
    fwrite(fileid,[' Name = "Body 1"' lb]);
    fwrite(fileid,[' Equation = 1' lb]);
    fwrite(fileid,[' Material = 1' lb]);
    fwrite(fileid,['Initial condition = 1' lb]);
    fwrite(fileid,['End' lb]);
    fwrite(fileid, lb);

    fwrite(fileid,['Body 2' lb]);
    fwrite(fileid,[' Target Bodies(1) = 2' lb]);
    fwrite(fileid,[' Name = "Body 2"' lb]);
    fwrite(fileid,[' Equation = 1' lb]);
    fwrite(fileid,[' Material = 2' lb]);
    fwrite(fileid,['End' lb]);
    fwrite(fileid, lb);

    fwrite(fileid,['Body 3' lb]);
    fwrite(fileid,[' Target Bodies(1) = 3' lb]);
    fwrite(fileid,[' Name = "Body 3"' lb]);
    fwrite(fileid,[' Equation = 1' lb]);
    fwrite(fileid,[' Material = 1' lb]);
    fwrite(fileid,['End' lb]);
    fwrite(fileid, lb);

%   Nonlinear Forced Vibration Analysis
%   ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    fwrite(fileid,['Solver 1' lb]);
    fwrite(fileid,['Equation = Nonlinear elasticity' lb]);
    fwrite(fileid,['Displace mesh = False' lb]);
    fwrite(fileid,[' Procedure = "ElasticSolve" "ElasticSolver"' lb]);
    fwrite(fileid,[' Element = p:2' lb]);    
    fwrite(fileid,['Variable = -dofs 3 Displacement' lb]);
    fwrite(fileid,[' Exec Solver = Always' lb]);
    fwrite(fileid,[' Stabilize = True' lb]);
    fwrite(fileid,[' Bubbles = False' lb]);
    fwrite(fileid,[' Lumped Mass Matrix = False' lb]);
    fwrite(fileid,[' Optimize Bandwidth = True' lb]);
    fwrite(fileid,['Displace mesh = False' lb]);
    fwrite(fileid,[' Steady State Convergence Tolerance = 1.0e-5' lb]);
    fwrite(fileid,[' Nonlinear System Convergence Tolerance = 1.0e-7' lb]);
    fwrite(fileid,[' Nonlinear System Max Iterations = 40' lb]);
    fwrite(fileid,[' Nonlinear System Newton After Iterations = 3' lb]);
    fwrite(fileid,[' Nonlinear System Newton After Tolerance = 1.0e-3' lb]);
    fwrite(fileid,[' Nonlinear System Relaxation Factor = 0.7' lb]);
    fwrite(fileid,[' Linear System Solver =  Direct !Iterative' lb]);
    fwrite(fileid,[' Linear System Iterative Method = BiCGStab' lb]);
    fwrite(fileid,[' Linear System Max Iterations = 10000' lb]);
    fwrite(fileid,['  Linear System Convergence Tolerance = 1.0e-7' lb]);
    fwrite(fileid,[' BiCGstabl polynomial degree = 2' lb]);
    fwrite(fileid,[' Linear System Preconditioning = ILU0' lb]);
    fwrite(fileid,[' Linear System ILUT Tolerance = 1.0e-3' lb]);
    fwrite(fileid,[' Linear System Abort Not Converged = False' lb]);
    fwrite(fileid,[' Linear System Residual Output = 10' lb]);
    fwrite(fileid,[' Linear System Precondition Recompute = 1' lb]);
    fwrite(fileid,['End' lb]);
    fwrite(fileid, lb);

    %Save some output files (for paraview %ect.)----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    fwrite(fileid,['Solver 2' lb]);
    fwrite(fileid,[' Exec Solver = After Simulation' lb]);
    fwrite(fileid,[' Equation = "result output"' lb]);
    fwrite(fileid,[' Procedure = "ResultOutputSolve" "ResultOutputSolver"' lb]);
    fwrite(fileid,[' Output File Name = "paraview"' lb]);
    fwrite(fileid,[' Binary Output = Logical True' lb]);
    fwrite(fileid,[' Single Precision = Logical True' lb]);
    %fwrite(fileid,[' Gid Format = Logical True' lb]);
    fwrite(fileid,[' Gmsh Format = Logical True' lb]);
    %fwrite(fileid,[' Vtk Format = Logical True' lb]);
    fwrite(fileid,[' Vtu Format = Logical True' lb]);
    %fwrite(fileid,[' Dx Format = Logical True' lb]);
    fwrite(fileid,[' Scalar Field 1 = String Potential' lb]);
    fwrite(fileid,['End' lb]);
    fwrite(fileid, lb);
    %*************
    fwrite(fileid,['Solver 3' lb]);
	fwrite(fileid,[' Exec Solver = After time step' lb]);
	fwrite(fileid,[' procedure = File "SaveData" "SaveScalars"' lb]);
	fwrite(fileid,[' Filename = "model_results.dat"' lb]);
	fwrite(fileid,[' Variable 1 = -dofs 3 Displacement' lb]);
% 	fwrite(fileid,[' Operator 1 = boundary sum' lb]);
	fwrite(fileid,['Save Coordinates(1,2,3) = Real 15.203 -0.01 9.7142' lb]) % Coordinate of points for results to be saved in the .dat file
	fwrite(fileid,['End' lb]);
	fwrite(fileid, lb);

    %Equation--------------------------------------------------------------------------------------------------------------------------------------------
    fwrite(fileid,['Equation 1' lb]);
    fwrite(fileid,[' Name = "Equation 1"' lb]);
    fwrite(fileid,[' Active Solvers(1) = 1' lb]);
    fwrite(fileid,['End' lb]);
    fwrite(fileid, lb);
    
    %Material properties for different parts of the structural member
    % Define Material properties 1-----------------------------------------------------------------------
    fwrite(fileid,['Material 1' lb]);
    fwrite(fileid,[' Name = "Material 1"' lb]);
    fwrite(fileid,[' Density = 2778' lb]);
    fwrite(fileid,[' Youngs modulus = 70.0e9' lb]);
    fwrite(fileid,[' Poisson ratio = 0.3' lb]);
    fwrite(fileid,[' Rayleigh Damping=Logical True' lb]);
    fwrite(fileid,[' Rayleigh Damping Beta=Real .01' lb]);
    fwrite(fileid,['End' lb]);
    fwrite(fileid, lb);

    % Define Material properties 2-----------------------------------------------------------------------
    fwrite(fileid,['Material 2' lb]);
    fwrite(fileid,[' Name = "Material 2"' lb]);
    fwrite(fileid,[' Density = 2778' lb]);
    fwrite(fileid,[' Youngs modulus = 51.03e9' lb]);
    fwrite(fileid,[' Poisson ratio = 0.3' lb]);
    fwrite(fileid,[' Rayleigh Damping=Logical True' lb]);
    fwrite(fileid,[' Rayleigh Damping Beta=Real .01' lb]);
    fwrite(fileid,['End' lb]);
    fwrite(fileid, lb);
 
   %boundary conditons------------------------------------------------------------------
 
   %boundary conditon 1 (for fixed condition of a cantilever beam)
    fwrite(fileid,['Boundary Condition 1' lb]);
    fwrite(fileid,['  Name = "BoundaryCondition 1"' lb]);
    fwrite(fileid,['  Target Boundaries(1) = 4' lb]);
    fwrite(fileid,['  Displacement 1 = 0' lb]);
    fwrite(fileid,['  Displacement 2 = 0' lb]);
    fwrite(fileid,['  Displacement 3 = 0' lb]);
    fwrite(fileid,['End' lb]);
    fwrite(fileid,lb);
    
    %Load input and boundary condition
    
    fwrite(fileid,['Boundary Condition 2' lb]);
    fwrite(fileid,['  Name = "BoundaryCondition 2"' lb]);
    fwrite(fileid,['  Target Boundaries(1) = 14' lb]);
    fwrite(fileid,['  Force 2 = Variable time' lb]);
    fwrite(fileid,['  Real MATC "750000*sin(41.88*tx(0))"' lb]); % Provide load in case of forced vibration
    fwrite(fileid,['End' lb]);
    fwrite(fileid,lb);

    %close file
    fclose(fileid);

    %Generate startinfo for elmersolver.exe=================================================================================
    %define path
    solverinfo_path =  [project_path filesep 'ELMERSOLVER_STARTINFO'];
 
   %Open file
    fileid = fopen(solverinfo_path,'w');

    %write file
    fwrite(fileid,['case.sif' lb]);
    fwrite(fileid,['1' lb])
    fwrite(fileid, lb)

    %close file
    fclose(fileid);

    %run simulation with elmersolver.exe==========================================================================================================
    
    disp("************************************************************************")
    simulation_command = [elmer_path filesep 'ElmerSolver.exe'];   
    system(simulation_command);

    %read results from fdisplacement_results.dat==============================================================================================
    
    fileid = fopen([project_path filesep 'model_results.dat'],'r');% name of the .dat file where intended coordinates results will be saved
    result_string = fgetl(fileid);
    fclose(fileid);
    
    %read results from paraview=======================================================================================
%     %calculate resistance and display results
%     paraview_command = ['D:\Elmer_Axial bar\New folder\case_t0001.vtu']
%     system(paraview_command);
%%% Section for post-processing 
    %read results from .dat=======================================================================================
    A=importdata('model_results.dat');
    dt=0.005;
    D=A(:,2);
    l=length(D);
    V=zeros(l,1);
    A=zeros(l,1);
    for i=2:l-1
        V(i,1)=(D(i+1,1)-D(i,1))/0.005; %% first order differentiation
        A(i,1)=(D(i+1,1)-(2*D(i,1))+D(i-1,1))/(0.005)^2; % 2nd order differentiation
    end
%%%****************************Develop Poincare' Map**********************************************************************    
Tspan=2.5; % Total time of analysis
t=(dt:0.005:Tspan)';
%lcv=3; % Provide least count of approximation for velocity; should be fixed based on the intended accuracy
multi_lcd=1e3; % multiplier to get integer number
lcd=1; % Provide least count of approximation for displacement;should be fixed based on the intended accuracy
st=1; %Provide assumed time for the initiation of stead-state
Et1=st/dt; %initiation point of steady-state
    pDisp=D(Et1:l,1); % 5001 
    pVel=V(Et1:l,1);
    pt=t(Et1:l,1);
    pn=zeros(size(pt,1),2);
    figure()
    set(gca,'fontWeight','bold','linewidth',1.5)
    grid on
    hold on
    itr = 1;

for j = round(length(pDisp(:,1))/2):length(pDisp(:,1))-1;
if(pDisp(j,1)>=pDisp((j-1),1))&&(pDisp(j,1)>=pDisp((j+1),1))
             pn(itr,1)=pDisp(j,1);
             pn(itr,2)=pVel(j,1);
             itr =itr+1;
            hold on;
        end
   end
    dis = pn(:,1).*multi_lcd;
    pn(:,1)=round(dis,lcd);
    [~,point] = unique(pn(:,1),'stable');
    poincare = pn(point,:);
    poincare( ~any(poincare,2), : ) = [];  
    poincare( :, ~any(poincare,1) ) = [];
    scatter(poincare(:,1).*1e-3,poincare(:,2),'xb','linewidth',1.1);
        hold off
    figure(3)
    plot(pDisp(:,1),pVel(:,1))
  

%%%%%%% Section to generate Poincare' map. More than one-points in th ePoincare' map either difines the presence of either transient part in the
%motion or sub or superharmonic frequency indicating a nonlinear system %%%%

  for i=Et1:length(t)
    subplot(2,3,1) 
    plot(t(i,1),D(i,1),'o','linewidth',2)
    xlabel('Time')
    ylabel('Displacement')
    hold on
    subplot(2,3,1) 
    plot(t(j:i),D(j:i),'Color','k','linewidth',2)
    xlabel('Time')
    ylabel('Displacement')
    title('Displacement time-History')   
       
    subplot(2,3,2) 
    plot(t(i,1),V(i,1),'o','linewidth',2)  
    xlabel('Time')
    ylabel('Velocity')
    hold on
    subplot(2,3,2)
    plot(t(j:i),V(j:i),'Color', 'k','linewidth',2)
    xlabel('Time')
    ylabel('Velocity')
    title('Velocity time-History')   
       
    subplot(2,3,3) 
    plot(t(i,1),A(i,1),'o','linewidth',2) 
    xlabel('Time')
    ylabel('Acceleration')
    hold on
    subplot(2,3,3) 
    plot(t(j:i),A(j:i),'Color', 'k','linewidth',2)
    xlabel('Time')
    ylabel('Acceleration')
    title('Acceleration time-History')   

    subplot(2,3,4) 
    plot(D(i,1),V(i,1),'o','linewidth',2)
    xlabel('Displacement')
    ylabel('Velocity')
    hold on
    subplot(2,3,4) 
    plot(D(j:i),V(j:i),'Color', 'k','linewidth',2)
    xlabel('Displacement')
    ylabel('Velocity')
    title('Phase-Space Diagram-2D')  

    subplot(2,3,5)
    plot3(t(i),D(i),V(i),'o','linewidth',2) %Poincare map
    xlabel('Time')
    ylabel('Displacement')
    zlabel('Velocity')
    hold on
    subplot(2,3,5) 
    plot3(t,D,V,'Color','k','linewidth',2)
    xlabel('Time')
    ylabel('Displacement')
    zlabel('Velocity')
    title('3D-Phase-Space Diagram and Poincare section') 
    
    subplot(2,3,6) 
    plot(poincare(:,1),poincare(:,2),'ro','linewidth',2) %Poincare map
    hold on
    xlabel('Displacement')
    ylabel('Velocity')
    title('Poincare map') 
    pause(0.05) 
    if i~=length(t)
        clf
    end
  end    
