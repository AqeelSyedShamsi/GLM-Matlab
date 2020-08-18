myfuns = moreVoxelBOLDstat;
pVal=[];
pValue=[];
DataModelallEffects=[];
Data_star=[];
stimDataX=[];
bVector=[];
designMatrix=[];
counfoundMatrix=[];
gaussianNoise=[];
sigmaMat=[];
modelDataX=[];

k=[0,1,2,3,4,5,6,7];
lVal=[0,1,2,3,4,5,6,7];
Mean =0; %Ground Truth mean
i = 300; %number of scans
j = 1; %number of voxels
% Creating a model of brain activity (BOLD)
res = 4;%Stimulus RT in seconds
%create a canonical BOLD response
TR=2;
HR=80/60;
BR=12/60;
alpha = 0.05;
loopItr = 10000;


gamma=[0.01:0.01:0.08];%;0.002:0.002:0.016;0.01:0.01:0.08];
GMean=0;
Gstd=[0.1];%,1,1.5];
pVal=[];
% for sni = 1:length(lVal)
%     l = lVal(sni);
for rhoi = 1:length(k)
    l = lVal(rhoi);
    T=[];
    betaMat=[];
for xx = 1:loopItr

    
    
    [X,G,H, Gnoise, resp, stim, hrf] = myfuns.glm_simulate_data(i,j,TR,HR,BR, gamma,res,GMean, k, l, rhoi);% X,G should be gernarate at once before the loop

    %Parameters 
    [b, beta, e, rss, r, sigma, inveG] = myfuns.glm_estimate(X,G);
    %betaMat(jx, xx)= b;      
    c=zeros(1,length(b));
    c(1)=1;
    c(2)=-1;
    betaMat=[betaMat b];
        
    t = myfuns.glm_inference (c, sigma, b, G,i,rss,inveG);
    %T (jx, xx)= t;
    %T=[T (t(1,:))'];
    T=[T t];
    

end
    [X_mod, X_star, modelDataAllEffects] = myfuns.glm_model(G,e,b,H,gamma,X, Gnoise);
    p = 1-tcdf(T,i-1);
    H1 = p < alpha;
    aq=sum(H1,2);
    aq1=sum(aq);
    disp(aq1/loopItr)
    pVal=[pVal (aq1/loopItr)];
    
    
    DataModelallEffects = [DataModelallEffects modelDataAllEffects];
    Data_star = [Data_star X_star];
    stimDataX= [stimDataX X];
    bVector=[bVector b];
    designMatrix=[designMatrix G];
    counfoundMatrix=[counfoundMatrix H];
    gaussianNoise=[gaussianNoise Gnoise];
    sigmaMat=[sigmaMat sigma];
    modelDataX=[modelDataX X_mod];
end
pValue=[pValue pVal];
% end