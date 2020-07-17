
myfuns = BOLDstat;


Mean =0; %Ground Truth mean

DataModelallEffects=[];
Data_star=[];
stimDataX=[];
bVector=[];
designMatrix=[];
counfoundMatrix=[];
gaussianNoise=[];

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

    
    
    [X,G,H, Gnoise, resp, stim, hrf] = myfuns.glm_simulate_data(i,j,TR,HR,BR, gamma,res,GMean, GStd);% X,G should be gernarate at once before the loop

    %Parameters 
    [b, beta, e, rss, r, sigma, inveG] = myfuns.glm_estimate(X,G);
    %betaMat(jx, xx)= b;      
    c=zeros(1,length(b));
    c(1)=1;
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
    disp(aq1/10000)
    pVal=[pVal (aq1/10000)];
    
    DataModelallEffects = [DataModelallEffects modelDataAllEffects];
    Data_star = [Data_star X_star];
    stimDataX= [stimDataX X];
    bVector=[bVector b];
    designMatrix=[designMatrix G];
    counfoundMatrix=[counfoundMatrix H];
    gaussianNoise=[gaussianNoise Gnoise];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Figure of BOLD Signal
figure('Name','BOLD Signal and Stiulus Response','NumberTitle','off')
subplot 311, plot(hrf),title('Response to a single stimulus')
subplot 312, plot([stim,resp]),title('A set of stimuli at randomized times')
subplot 313, plot(resp),title('Response to the whole set of stimuli')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Data Structre graphs separate for different variances
figure('Name','Guassian Noise for different vaiances','NumberTitle','off')
subplot(3,2,1), plot(gaussianNoise(:,1)); title('Gaussian noise with std of 0.1')
subplot(3,2,2),hist(gaussianNoise(:,1),20);title('Gaussian noise histogram with std of 0.1');
subplot(3,2,3), plot(gaussianNoise(:,2)); title('Gaussian noise with std of 1')
subplot(3,2,4),hist(gaussianNoise(:,2),20);title('Gaussian noise histogram with std of 1');
subplot(3,2,5), plot(gaussianNoise(:,3)); title('Gaussian noise with std of 1.5')
subplot(3,2,6),hist(gaussianNoise(:,3),20);title('Gaussian noise  histogram with std of 1.5');

figure('Name','Counfound Data for different vaiances','NumberTitle','off')
subplot(3,2,1), plot(counfoundMatrix(:,[1:8])); title('Confound Data with std of 0.1')
subplot(3,2,2),hist(counfoundMatrix(:,[1:8]),20);title('Confound Data histogram with std of 0.1');
subplot(3,2,3), plot(counfoundMatrix(:,[9:16])); title('Confound Data with std of 1')
subplot(3,2,4),hist(counfoundMatrix(:,[9:16]),20);title('Confound Data histogram with std of 1');
subplot(3,2,5), plot(counfoundMatrix(:,[17:24])); title('Confound Data with std of 1.5')
subplot(3,2,6),hist(counfoundMatrix(:,[17:24]),20);title('Confound Data  histogram with std of 1.5');

figure('Name','Design Matrix for different vaiances','NumberTitle','off')
subplot(3,2,1), plot(designMatrix(:,[1:10])); title('Design Matrix with std of 0.1')
subplot(3,2,2),imagesc(designMatrix(:,[1:10]));title('Design Matrix heatmap with std of 0.1');
subplot(3,2,3), plot(designMatrix(:,[11:20])); title('Design Matrix with std of 1')
subplot(3,2,4),imagesc(designMatrix(:,[11:20]));title('Design Matrix heatmap with std of 1');
subplot(3,2,5), plot(designMatrix(:,[21:30])); title('Design Matrix with std of 1.5')
subplot(3,2,6),imagesc(designMatrix(:,[21:30]));title('Design Matrix  heatmap with std of 1.5');

figure('Name','Simulated Data for different vaiances','NumberTitle','off')
subplot(3,2,1), plot(stimDataX(:,1)); title('Simulated Data with std of 0.1')
subplot(3,2,2),hist(stimDataX(:,1),20);title('Simulated Data histogram with std of 0.1');
subplot(3,2,3), plot(stimDataX(:,2)); title('Simulated Data with std of 1')
subplot(3,2,4),hist(stimDataX(:,2),20);title('Simulated Data histogram with std of 1');
subplot(3,2,5), plot(stimDataX(:,3)); title('Simulated Data with std of 1.5')
subplot(3,2,6),hist(stimDataX(:,3),20);title('Simulated Data  histogram with std of 1.5');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Data Structre graphs together for different variances
figure('Name','Data Structure','NumberTitle','off')
% subplot(421), plot(Gnoise); title('Gaussian noise')
% subplot(422),hist(Gnoise,20);title('Gaussian noise histogram');
% subplot(423) ,plot(H); title('Counfound Matrix')
% subplot(424),hist(H,20);title('Counfound Matrix histogram');
% subplot(425) ,plot(G); title('Design Matrix')
% subplot(426),imagesc(G);title('Design Matrix histogram');
% subplot(427) ,plot(X); title('Data Orignal')
% subplot(428),hist(X,20);title('Data Matrix histogram');
subplot(12,2,1), plot(gaussianNoise(:,1)); title('Gaussian noise with std of 0.1')
subplot(12,2,2),hist(gaussianNoise(:,1),20);title('Gaussian noise histogram with std of 0.1');
subplot(12,2,3), plot(gaussianNoise(:,2)); title('Gaussian noise with std of 1')
subplot(12,2,4),hist(gaussianNoise(:,2),20);title('Gaussian noise histogram with std of 1');
subplot(12,2,5), plot(gaussianNoise(:,3)); title('Gaussian noise with std of 1.5')
subplot(12,2,6),hist(gaussianNoise(:,3),20);title('Gaussian noise  histogram with std of 1.5');

subplot(12,2,7), plot(counfoundMatrix(:,[1:8])); title('Confound Data with std of 0.1')
subplot(12,2,8),hist(counfoundMatrix(:,[1:8]),20);title('Confound Data histogram with std of 0.1');
subplot(12,2,9), plot(counfoundMatrix(:,[9:16])); title('Confound Data with std of 1')
subplot(12,2,10),hist(counfoundMatrix(:,[9:16]),20);title('Confound Data histogram with std of 1');
subplot(12,2,11), plot(counfoundMatrix(:,[17:24])); title('Confound Data with std of 1.5')
subplot(12,2,12),hist(counfoundMatrix(:,[17:24]),20);title('Confound Data  histogram with std of 1.5');

subplot(12,2,13), plot(designMatrix(:,[1:10])); title('Design Matrix with std of 0.1')
subplot(12,2,14),imagesc(designMatrix(:,[1:10]));title('Design Matrix heatmap with std of 0.1');
subplot(12,2,15), plot(designMatrix(:,[11:20])); title('Design Matrix with std of 1')
subplot(12,2,16),imagesc(designMatrix(:,[11:20]));title('Design Matrix heatmap with std of 1');
subplot(12,2,17), plot(designMatrix(:,[21:30])); title('Design Matrix with std of 1.5')
subplot(12,2,18),imagesc(designMatrix(:,[21:30]));title('Design Matrix  heatmap with std of 1.5');

subplot(12,2,19), plot(stimDataX(:,1)); title('Simulated Data with std of 0.1')
subplot(12,2,20),hist(stimDataX(:,1),20);title('Simulated Data histogram with std of 0.1');
subplot(12,2,21), plot(stimDataX(:,2)); title('Simulated Data with std of 1')
subplot(12,2,22),hist(stimDataX(:,2),20);title('Simulated Data histogram with std of 1');
subplot(12,2,23), plot(stimDataX(:,3)); title('Simulated Data with std of 1.5')
subplot(12,2,24),hist(stimDataX(:,3),20);title('Simulated Data  histogram with std of 1.5');



figure('Name','T Statistics','NumberTitle','off')
bar(pVal);
xlabel('Gamma');
ylabel('alpha'); 




[X_mod, X_star, modelDataAllEffects] = myfuns.glm_model(G,e,b,H,gamma,X, Gnoise);

% figure('Name','Data with all Effects','NumberTitle','off')
% subplot 311, plot(X),title('Orignal Data')
% subplot 312, plot(modelDataAllEffects),title('model data with all effects')
% subplot 313, plot([X,modelDataAllEffects]),title('fitted model'); legend('Orignal Data','Modeled Data');

figure('Name','Data with all Effects','NumberTitle','off')
subplot 621, plot(stimDataX(:,1)),title('Simulated Data (low gamma and std of 0.1)')
subplot 622, plot(DataModelallEffects(:,1)),title('model data with effects of interest')
subplot 612, plot([stimDataX(:,1),DataModelallEffects(:,1)]),title('fitted model'); legend('Simulated Data','Data with effects of interest');
subplot 625, plot(stimDataX(:,2)),title('Simulated Data (medium gamma and std of 1')
subplot 626, plot(DataModelallEffects(:,2)),title('model data with effects of interest')
subplot 614, plot([stimDataX(:,2),DataModelallEffects(:,2)]),title('fitted model'); legend('Simulated Data','Data with effects of interest');
subplot 629, plot(stimDataX(:,3)),title('Simulated Data (high gamma and std of 1.5')
subplot (6,2,10), plot(DataModelallEffects(:,3)),title('model data with effects of interest')
subplot 616, plot([stimDataX(:,3),DataModelallEffects(:,3)]),title('fitted model'); legend('Simulated Data','Data with effects of interest');
%subplot 313, plot(X,modelDataAllEffects),title('Linear Regression')

figure('Name','Adjusted Data and model using all effects','NumberTitle','off')
subplot 621, plot(Data_star(:,1)),title('Adjusted Data (low gamma and std of 0.1)')
subplot 622, plot(DataModelallEffects(:,1)),title('model data with effects of interest')
subplot 612, plot([Data_star(:,1),DataModelallEffects(:,1)]),title('fitted model'); legend('Adjusted Data','Data with effects of interest');
subplot 625, plot(Data_star(:,2)),title('Adjusted Data (medium gamma and std of 1)')
subplot 626, plot(DataModelallEffects(:,2)),title('model data with effects of interest')
subplot 614, plot([Data_star(:,2),DataModelallEffects(:,2)]),title('fitted model'); legend('Adjusted Data','Data with effects of interest');
subplot 629, plot(Data_star(:,3)),title('Adjusted Data (high gamma and std of 1.5)')
subplot (6,2,10), plot(DataModelallEffects(:,3)),title('model data with effects of interest')
subplot 616, plot([Data_star(:,3),DataModelallEffects(:,3)]),title('fitted model'); legend('Adjusted Data','Data with effects of interest');


figure('Name','Beta Plot with different variances','NumberTitle','off')
plot(bVector),title('Beta Plot'); legend('bVector 0.1','bVector 1','bVector 1.5')

% figure('Name','Adjusted Data and model using all effects','NumberTitle','off')
% subplot 311, plot(X_star),title('Adjusted Data')
% subplot 312, plot(X_mod),title('model data with effects of interest')
% subplot 313, plot([X_star,X_mod]),title('fitted model'); legend('Adjusted Data','Data with effects of interest');


