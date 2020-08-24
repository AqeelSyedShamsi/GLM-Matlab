myfuns = FTestStat;

i = 300;                            %number of scans
j = 1;                              %number of voxels
% Creating a model of brain activity (BOLD)
res = 4;                            %Stimulus RT in seconds
%create a canonical BOLD response
TR=2;

alpha = 0.05;
loopItr = 10000;
F=[];
pVal=[];

for xx = 1:loopItr

    
    [X_FM,X_RM, G_FM,G_RM, resp] = myfuns.glm_simulate_data(i,res);
    %Parameters
    [b_FM, b_RM, r_FM, r_RM,e_FM,e_RM] =myfuns.glm_estimate(G_FM,G_RM, X_FM, X_RM);
   
    [X_PFM, X_PRM,SSE_FM, SSE_RM] =myfuns.glm_model(G_FM,G_RM, b_FM, b_RM, e_FM, e_RM);  
    F_star = myfuns.glm_inference (SSE_RM, SSE_FM, r_RM, r_FM); %(c, beta, sigma, alpha)
    F=[F F_star];
    

end
p = 1-fcdf(F,r_RM,r_FM);
H1 = p < alpha;
aq=sum(H1,2);
aq1=sum(aq);
disp(aq1/loopItr)