%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GLM for BOLD simulation with two different conditions%

myfuns = voxelBOLDstat22_07;

pVal_FP=[];                         %false positive 
pVal=[];                            %p value             
CmulB=[];                           %c*b
DataModelallEffects=[];              
Data_star=[];                       %adjusted data
stimDataX=[];
bVector=[];
designMatrix=[];
counfoundMatrix=[];
gaussianNoise=[];
sigmaMat=[];
modelDataX=[];
k=[0,1,2,3,4,5,6,7];                %rhoi function 
lVal=[0,1,2,3,4,5,6,7];             %sn function
Mean =0;                            %Ground Truth mean
i = 300;                            %number of scans
j = 1;                              %number of voxels
% Creating a model of brain activity (BOLD)
res = 4;                            %Stimulus RT in seconds
%create a canonical BOLD response
TR=2;
HR=80/60;                           %Heart rate per minute
BR=12/60;                           %breath rate per minute
alpha = 0.05;
loopItr = 10000;


gamma=[0.01:0.01:0.08];%;0.002:0.002:0.016;0.01:0.01:0.08];
GMean=0;                            %mean for gaussian noise
%Gstd=[0.1];%,1,1.5];               %std of gaussian noise
pVal=[];                            %vector for p values
TMat=[];
for rhoi = 1:length(k)              %loop for different k and l values
    l = lVal(rhoi);
    T=[];
    betaMat=[];
for xx = 1:loopItr

    
    
    [X,G,H, Gnoise, resp, hrf] = myfuns.glm_simulate_data(i,TR,HR,BR, gamma,res,GMean, k, l, rhoi);
    %Parameters 
    [b, beta, e, rss, r, sigma, inveG] = myfuns.glm_estimate(X,G);
    %betaMat(jx, xx)= b;      
    c=zeros(1,length(b));
    c(1)=1;
    c(2)=-1;
    betaMat=[betaMat b];
        
    t = myfuns.glm_inference (c, b);
    
    T=[T t];
    

end
    [X_mod, X_star, modelDataAllEffects] = myfuns.glm_model(G,b,H,X, Gnoise);
    p = 1-tcdf(T,i-1);
    
    H1 = p < alpha;
    H1_FP=p>alpha;
    aq=sum(H1,2);
    aq1=sum(aq);
    disp(aq1/loopItr)
    pVal=[pVal (aq1/loopItr)];
    pVal_FP=[pVal_FP ((sum(aq))/loopItr) ];
    
    
    
    DataModelallEffects = [DataModelallEffects modelDataAllEffects];
    Data_star = [Data_star X_star];
    stimDataX= [stimDataX X];
    bVector=[bVector b];
    designMatrix=[designMatrix G];
    counfoundMatrix=[counfoundMatrix H];
    gaussianNoise=[gaussianNoise Gnoise];
    sigmaMat=[sigmaMat sigma];
    modelDataX=[modelDataX X_mod];
    CmulB=[CmulB (c*b)];
end

figure('Name','T Statistics','NumberTitle','off')
bar(pVal);
xlabel('Gamma');
ylabel('alpha'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Figure of BOLD Signal
figure('Name','BOLD Signal and Stiulus Response','NumberTitle','off')
subplot 311, plot(hrf),title('Response to a single stimulus')
subplot 312, plot(resp),title('Response to the whole set of stimuli')

figure('Name','Data with all Effects','NumberTitle','off')
subplot 621, plot(stimDataX(:,1)),title('Simulated Data)')
subplot 622, plot(DataModelallEffects(:,1)),title('model data with effects of interest')
subplot 612, plot([stimDataX(:,1),DataModelallEffects(:,1)]),title('fitted model'); legend('Simulated Data','Data with effects of interest');
subplot 625, plot(stimDataX(:,2)),title('Simulated Data')
subplot 626, plot(DataModelallEffects(:,2)),title('model data with effects of interest')
subplot 614, plot([stimDataX(:,2),DataModelallEffects(:,2)]),title('fitted model'); legend('Simulated Data','Data with effects of interest');
subplot 629, plot(stimDataX(:,3)),title('Simulated Data')
subplot (6,2,10), plot(DataModelallEffects(:,3)),title('model data with effects of interest')
subplot 616, plot([stimDataX(:,3),DataModelallEffects(:,3)]),title('fitted model'); legend('Simulated Data','Data with effects of interest');
%subplot 313, plot(X,modelDataAllEffects),title('Linear Regression')

figure('Name','Adjusted Data and model using all effects','NumberTitle','off')
subplot 621, plot(Data_star(:,6)),title('Adjusted Data')
subplot 622, plot(DataModelallEffects(:,6)),title('model data with effects of interest')
subplot 612, plot([Data_star(:,6),DataModelallEffects(:,6)]),title('fitted model'); legend('Adjusted Data','Data with effects of interest');
subplot 625, plot(Data_star(:,7)),title('Adjusted Data')
subplot 626, plot(DataModelallEffects(:,7)),title('model data with effects of interest')
subplot 614, plot([Data_star(:,7),DataModelallEffects(:,7)]),title('fitted model'); legend('Adjusted Data','Data with effects of interest');
subplot 629, plot(Data_star(:,8)),title('Adjusted Data')
subplot (6,2,10), plot(DataModelallEffects(:,8)),title('model data with effects of interest')
subplot 616, plot([Data_star(:,8),DataModelallEffects(:,8)]),title('fitted model'); legend('Adjusted Data','Data with effects of interest');


figure('Name','Beta Plot with different k and l values','NumberTitle','off')
plot(bVector),title('Beta Plot'); % for different k and l values
grid on;

figure('Name','Beta contrast c’*beta','NumberTitle','off')
plot(CmulB),title('Beta Plot');  % c*b for k and l values
grid on;

figure('Name','False positive rate','NumberTitle','off')
plot(pVal_FP),title('False Positive Rate');  % c*b for k and l values
grid on;


