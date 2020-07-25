myfuns = moreVoxelBOLDstat;
Mean =0; %Ground Truth mean

DataModelallEffects=[];
Data_star=[];
stimDataX=[];
bVector=[];
designMatrix=[];
counfoundMatrix=[];
gaussianNoise=[];
sigmaMat=[];
modelDataX=[];

i = 300; %number of scans
j = 1; %number of voxels

% Creating a model of brain activity (BOLD)
res = 2.4;%Stimulus RT in seconds
%duration = 20; % seconds
% create a canonical BOLD response
TR=2;
HR=80/60;
BR=12/60;
alpha = 0.05;
loopItr = 10000;
% betaMat=zeros (j,loopItr); %beta matrix
%T=zeros (j,loopItr);       %T test matrix
T=[];
betaMat=[];
%design matrix
% c=zeros(1,j);
% c(1)=1;
%c=[1,0,0,0,0,0,0];
%c=c';
%creating gaussian signal and Design matrix
gama=[0.001:0.001:0.008;0.002:0.002:0.016;0.01:0.01:0.08];
GMean=0;
Gstd=[0.1,1,1.5];
pVal=[];

for gamaX = 1:length(gama(:,2))
    gamma = gama(gamaX,:);
    
    GStd  = Gstd(gamaX);
    for xx = 1:loopItr

    
    
        [X,G,H, Gnoise, resp, hrf] = myfuns.glm_simulate_data(i,j,TR,HR,BR, gamma,res,GMean, GStd);% X,G should be gernarate at once before the loop

        %Parameters 
        [b, beta, e, rss, r, sigma, inveG] = myfuns.glm_estimate(X,G);
        %betaMat(jx, xx)= b;      
        c=zeros(1,length(b));
        c(1)=1;
        c(2)=-1;
        betaMat=[betaMat b];
        t = myfuns.glm_inference (c, b);
        %T (jx, xx)= t;
        %T=[T (t(1,:))'];
        T=[T t'];
        
   
    end
    p = 1-tcdf(T,i-1);
    H1 = T < alpha;
    aq=sum(H1,2);
    aq1=sum(aq);
    disp(aq1/30000)
    pVal=[pVal (aq1/30000)];
end