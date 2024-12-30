%% PSO Satellite Orbit Generator
% This code simply generates desired random number of satellites, Orbits and stations
% and selects best orbits using PSO algorithm. 
% 'samp' is number of orbits to be generate
% 'nf' is number of best orbits


%% Create a Satellite Scenario
clear;
sc = satelliteScenario;
load("timetableSatelliteTrajectory.mat","positionTT","velocityTT");
% Sat Params
startTime = datetime(2020, 5, 1, 11, 36, 0);
stopTime = startTime + days(1);
sampleTime = 60;
sc = satelliteScenario(startTime, stopTime, sampleTime);

%% Creating satellites and ground stations
% Number of samples to generate
samp=50;
%
% latitudes
a = 1;
b = 89;
r = (b-a).*rand(samp,1) + a;
r=round(r);
% longitudes
aa = -100;
bb = 1000;
rr = (bb-aa).*rand(samp,1) + aa;
rr=round(rr);
% semiMajorAxis 
aaa = 9000000;
bbb = 90000000;
rrr = (bbb-aaa).*rand(samp,1) + aaa;
rrr=round(rrr);
% eccentricity
aaaa = 0.1;
bbbb = 0.9;
rrrr = (bbbb-aaaa).*rand(samp,1) + aaaa;
% inclination
aaaaa = 5;
bbbbb = 360;
rrrrr = (bbbbb-aaaaa).*rand(samp,1) + aaaaa;
rrrrr=round(rrrrr);
% Final matrix to select best orbits
finalr=[r rr rrr rrrr rrrrr]';

%% PSO 
% Data Preparation
x=finalr(:,1:end-1)';
t=finalr(:,end)';
data.x=x;
data.t=t;
data.nx=size(x,1);
data.nt=size(t,1);
data.nSample=size(x,2);

%% Number of Desired PSO Orbits
nf=7;

%% Cost Function
CostFunction=@(u) FeatureSelectionCost(u,nf,data);
% Number of Decision Variables
nVar=data.nx;
% Size of Decision Variables Matrix
VarSize=[1 nVar];
% Lower Bound of Variables
VarMin=0;
% Upper Bound of Variables
VarMax=1;

%% PSO Parameters
% Maximum Number of Iterations
MaxIt=6;
% Population Size (Swarm Size)
nPop=3;
% Constriction Coefficients
phi1=2.05;
phi2=2.05;
phi=phi1+phi2;
chi=2/(phi-2+sqrt(phi^2-4*phi));
% Inertia Weight
w=chi;
% Inertia Weight Damping Ratio
wdamp=1;
% Personal Learning Coefficient
c1=chi*phi1;
% Global Learning Coefficient
c2=chi*phi2;
% Velocity Limits
VelMax=0.1*(VarMax-VarMin);
VelMin=-VelMax;

%% Basics
empty_particle.Position=[];
empty_particle.Cost=[];
empty_particle.Out=[];
empty_particle.Velocity=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];
empty_particle.Best.Out=[];
particle=repmat(empty_particle,nPop,1);
BestSol.Cost=inf;
for i=1:nPop
% Begin Position
particle(i).Position=unifrnd(VarMin,VarMax,VarSize);
% Begin Velocity
particle(i).Velocity=zeros(VarSize);
% Evaluation
[particle(i).Cost, particle(i).Out]=CostFunction(particle(i).Position);
% Update Personal Best
particle(i).Best.Position=particle(i).Position;
particle(i).Best.Cost=particle(i).Cost;
particle(i).Best.Out=particle(i).Out;
% Update Global Best
if particle(i).Best.Cost<BestSol.Cost
BestSol=particle(i).Best;
end
end
%
BestCost=zeros(MaxIt,1);

%% PSO Body Part

for it=1:MaxIt
for i=1:nPop
% Update Velocity
particle(i).Velocity = w*particle(i).Velocity ...
+c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
+c2*rand(VarSize).*(BestSol.Position-particle(i).Position);
% Apply Velocity Limits
particle(i).Velocity = max(particle(i).Velocity,VelMin);
particle(i).Velocity = min(particle(i).Velocity,VelMax);
% Update Position
particle(i).Position = particle(i).Position + particle(i).Velocity;
% Velocity Mirror Effect
IsOutside=(particle(i).Position<VarMin | particle(i).Position>VarMax);
particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
% Apply Position Limits
particle(i).Position = max(particle(i).Position,VarMin);
particle(i).Position = min(particle(i).Position,VarMax);
% Evaluation
[particle(i).Cost, particle(i).Out] = CostFunction(particle(i).Position);
% Update Personal Best
if particle(i).Cost<particle(i).Best.Cost
particle(i).Best.Position=particle(i).Position;
particle(i).Best.Cost=particle(i).Cost;
particle(i).Best.Out=particle(i).Out;
% Update Global Best
if particle(i).Best.Cost<BestSol.Cost
BestSol=particle(i).Best;
end
end
end
%
BestCost(it)=BestSol.Cost;
%
disp(['In Iteration ' num2str(it) ': PSO Fittest Value Is : ' num2str(BestCost(it))]);
w=w*wdamp;
end

%% Final Orbits
% Extracting Data
RealData=data.x';
FinalFeaturesInd=BestSol.Out.S;
% Sort
PSO_Orbits=RealData(:,FinalFeaturesInd);
%
% Ground Stations
for i=1:nf
gs(i) = groundStation(sc, PSO_Orbits(1,i), PSO_Orbits(2,i));
end;
% Add satellites using Keplerian elements.
rightAscensionOfAscendingNode = 0; 
argumentOfPeriapsis = 0; 
trueAnomaly = 0; 
for i=1:nf
satSGP4(i)= satellite(sc, PSO_Orbits(3,i), PSO_Orbits(4,i), PSO_Orbits(5,i), ...
        rightAscensionOfAscendingNode, argumentOfPeriapsis, trueAnomaly);
end;
% Result
play(sc);






