myfuns = stat;
mean =0;
i = 300; %number of repetations
j = 1; %1000; %number of voxels

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
gama=[0.1:0.1:0.8;1:1:8;10:10:80];
pVal=[];

for gamaX = 1:length(gama(:,2))
    gamma = gama(gamaX,:);
for xx = 1:loopItr

    
    
    [X,G,H] = myfuns.glm_simulate_data(i,j,TR,HR,BR, gamma);% X,G should be gernarate at once before the loop

    %Parameters 
    [b, beta, e, rss, r, sigma, inveG] = myfuns.glm_estimate(X,G);
    %betaMat(jx, xx)= b;      
    c=zeros(1,length(b));
    c(1)=1;
    betaMat=[betaMat b];
        
    t = myfuns.glm_inference (c, sigma, b, G,i,rss,inveG);
    %T (jx, xx)= t;
    T=[T (t(1,:))'];
    
    

end
    p = 1-tcdf(T,i-1);
    H1 = p < alpha;
    aq=sum(H1,2);
    aq1=sum(aq);
    disp(aq1/10000)
    pVal=[pVal (aq1/10000)];
    
end
bar(pVal);
xlabel('Gamma');
ylabel('alpha'); 
