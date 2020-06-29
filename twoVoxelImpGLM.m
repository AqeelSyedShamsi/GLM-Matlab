myfuns = withoutConfoundstat;
mean =0;
i = 300; %number of repetations
j = 5; %1000; %number of voxels


alpha = 0.05;
loopItr = 10000;
% betaMat=zeros (j,loopItr); %beta matrix
% T=zeros (j,loopItr);       %T test matrix
T=[];
betaMat=[];
%design matrix
c=zeros(1,j);
c(1)=1;
% c=1;
%creating gaussian signal and Design matrix
    

    
for xx = 1:loopItr

[X,G] = myfuns.glm_simulate_data(i,j);% X,G should be gernarate at once before the loop

    

%Parameters 
[b, beta, e, rss, r, sigma, inveG] = myfuns.glm_estimate(X,G);
%G = myfuns.glm_simulate_data(i,j);
%betaMat(jx, xx)= b;      

betaMat=[betaMat (b(1,:))'];
        
t = myfuns.glm_inference (c, sigma, b, G,i,rss,inveG);
%T (jx, xx)= t;
T=[T (t(1,:))'];
    

end


p = 1-tcdf(T,i-1);
H1 = p < alpha;
aq=sum(H1,2);
aq1=sum(aq);
disp(aq1/10000)



% p = 1-tcdf(T(4,:),i-1);
% H1 = p < alpha;
% aq=sum(H1,2);
% aq1=sum(aq);
% disp(aq1/10000)
% 

