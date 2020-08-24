function func = FTestStat
   
    func.glm_simulate_data=@glm_simulate_data;
    func.glm_estimate=@glm_estimate;
    func.glm_inference=@glm_inference;
    func.glm_model=@glm_model;

end
function [X_FM,X_RM, G_FM,G_RM, resp] = glm_simulate_data(i,res)
    G_FM=ones(i,1); %design matrix
    G_RM=ones(i,1); %design matrix
    
    hrf = spm_hrf(1/(res));
    hrfCod2 = spm_hrf(1/(res+(4/10)));%k(rhoi)+1*(hrf);
   
    resp=[];
    resp_RM=[];

    hrf1=hrf(1:30);   %reshaping of the signal
    hrf2=hrfCod2(1:30);%reshaping of the signal
    desCondA= zeros(300,1);
    desCondB= zeros(300,1);
    for respi = 1:i/60
        resp=vertcat(resp, hrf1);
        resp=vertcat(resp,hrf2);
        %resp_RM=vertcat(resp_RM, hrf1);
        
        
    end
    for respii = 1:i/30
        resp_RM=vertcat(resp_RM, hrf1);    
    end
    
    desCondA(1:30)=hrf1;
    desCondA(61:90)=hrf1;
    desCondA(121:150)=hrf1;
    desCondA(181:210)=hrf1;
    desCondA(241:270)=hrf1;
    
    desCondB(31:60)=hrf2;
    desCondB(91:120)=hrf2;
    desCondB(151:180)=hrf2;
    desCondB(211:240)=hrf2;
    desCondB(271:300)=hrf2;
    
    G_FM= [G_FM desCondA desCondB];
    X_FM= resp;
    X_RM= resp_RM;
    G_RM= [G_RM resp_RM];
    
end
function [b_FM, b_RM, r_FM, r_RM,e_FM,e_RM] =glm_estimate(G_FM,G_RM, X_FM, X_RM)
    inveG_FM=inv(G_FM'*G_FM);
    inveG_RM=inv(G_RM'*G_RM);


    
    b_FM= (inveG_FM*G_FM')*X_FM;
    b_RM= (inveG_RM*G_RM')*X_RM;
    %beta = mean(b);
    e_FM = X_FM-G_FM*b_FM;
    e_RM = X_RM-G_RM*b_RM;

    rss_FM = e_FM'*e_FM;
    rss_RM = e_RM'*e_RM;
    r_FM=size(X_FM,1)-rank(G_FM);
    r_RM=size(X_RM,1)-rank(G_RM);
    %sigma= sqrt(rss/r);
    
end

function [X_PFM, X_PRM,SSE_FM, SSE_RM] =glm_model(G_FM,G_RM, b_FM, b_RM, e_FM, e_RM)
    
    
    
    X_PFM= G_FM*b_FM;% + e
    X_PRM= G_RM*b_RM;% + e
    e_PFM = X_PFM-G_FM*b_FM;
    e_PRM = X_PRM-G_RM*b_RM;
    
    
    SSE_FM=sum((e_FM-e_PFM).^2);
    SSE_RM=sum((e_RM-e_PRM).^2);
    
    
    
end

function [F_star] = glm_inference (SSE_RM, SSE_FM, r_RM, r_FM) %(c, beta, sigma, alpha)

    F_star = (((SSE_RM - SSE_FM)/(r_RM - r_FM)))/(SSE_FM/r_FM);
end
