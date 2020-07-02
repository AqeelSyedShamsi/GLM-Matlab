myfuns = withoutConfoundstat;
mean =0;
i = 300; %number of repetations
j = 1; %1000; %number of voxels


alpha = 0.05;
loopItr = 10000;
betaMat=zeros (j,loopItr); %beta matrix
T=zeros (j,loopItr);       %T test matrix
%design matrix
c=zeros(1,j);
c(1)=1;
%creating gaussian signal and Design matrix
    

for xx = 1:loopItr

    [X,G] = myfuns.glm_simulate_data(i,j);% X,G should be gernarate at once before the loop

    

    %Parameters 
    [b, beta, e, rss, r, sigma, inveG] = myfuns.glm_estimate(X,G);
    %G = myfuns.glm_simulate_data(i,j);
    
    t = myfuns.glm_inference (c, sigma, b, G,i,rss,inveG);
    T (j, xx)= t;
    

end

p = 1-tcdf(T,i-1);
H1 = p < alpha;
disp(sum(H1)/10000)
