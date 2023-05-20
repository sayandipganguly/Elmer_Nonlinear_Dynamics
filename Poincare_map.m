%{
Script for developing Poincare map from a signal.

INPUT:
Steady-state Displacement, Velocity and acceleration time-history 

OUTPUT:
   Poincare Points and Poincare map

PROJECT: Research in Structural Health Monitoring - Indian Institute of
Technology Patna

AUTHOR: Sayandip Ganguly
%}
lcv=5; % Least count of approximation for velocity; should be fixed based on the intended accuracy
lcd=5; % Least count of approximation for displacement;should be fixed based on the intended accuracy
t=xlsread('result_node_11.xlsx','A2:A7148'); % Time e.g. 0:0.005:20
A=xlsread('result_node_11.xlsx','D2:D7148'); % Measured Acceleration of specified node of the member
V=xlsread('result_node_11.xlsx','C2:C7148');% Measured Velocity of specified node of the member
D=xlsread('result_node_11.xlsx','B2:B7148');% Measured Displacement of specified node of the member
Acc=0.5*max(A); % Point from acceleration time-history which will be checcked for repetition of specified node of the member
L2=length(A);
point=Acc*ones(L2,1)';
[t01,A01] = intersections(t,A,t,point,1); % Returns intersection points x02=Time, y02=intended displacement where intersection is required
figure(1)
plot(t,A,t,point,t01,A01,'ok');
L12=length(t01);
k1=1;
k2=1;
for i2=1:((L12-1)/2)
    Pvel2(i2)=round(interp1(t,V,t01(k1)),lcv); % Returns velocity at intersection times
    k1=k1+2;
end
for j2=1:((L12-1)/2)
    Pdis2(j2)=round(interp1(t,D,t01(k2)),lcd); % Returns Dsiplacement at intersection times
    k2=k2+2;
end
figure(1) 
plot(Pdis2,Pvel2,'r.') %Poincare map
figure(2)
plot(t,A,t,point,t01,A01,'ok');