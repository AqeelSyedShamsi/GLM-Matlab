function func = voxelBOLDstat22_07
   
    func.glm_simulate_data=@glm_simulate_data;
    func.glm_estimate=@glm_estimate;
    func.glm_inference=@glm_inference;
    func.glm_model=@glm_model;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Simulation Function%%%%%%%%%%%%%%%%
function [X, G, H, Gnoise, resp, hrf] = glm_simulate_data(i,TR,HR,BR, gamma,res,GMean, k, l, rhoi)
    
    G=ones(i,1); %design matrix
    Gnoise = GMean + (0.5/(l+1)).*randn(i,1);%Gaussian Noise with mean 0 and standard deviation wrt l.

       
    hrf = spm_hrf(1/(res));%condition A BOLD
    hrfCod2 = spm_hrf(1/(res+(k(rhoi)/10)));%condition B BOLD
   
    resp=[];

    hrf1=hrf(1:30);   %reshaping of the signal
    hrf2=hrfCod2(1:30);%reshaping of the signal
    desCondA= zeros(300,1);
    desCondB= zeros(300,1);
    for respi = 1:i/60
        resp=vertcat(resp, hrf1);
        resp=vertcat(resp,hrf2);
        
        
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

    
    
    %n1=1:1:i;
    %n=linspace (0,0.001,300) ;
    n=0.01:0.01:3;
    t=n*TR;
    hrFreq=sin(2*pi*HR*t); %Heart signal
    brFreq=sin(2*pi*BR*t); %Breathing signal

    
    inputdata =rand(i,6); %random data for confounds
    H= (random('normal',0,(0.5/(l+1)),size(inputdata))); %RANDOM data
    H =[H hrFreq'];
    H =[H brFreq'];
    
    G= [G desCondA desCondB H]; %Design matrix with BOLD and confounding effects
    X= resp+ Gnoise + (H*gamma'); %Simulated Data
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%Estimation functions%%%%%%%%%%
function [b, beta, e, rss, r, sigma, inveG] =glm_estimate(X,G)
    inveG=inv(G'*G);
    %b = pinv(G)*X;
    
    b= (inveG*G')*X; %beta parameters
   
    beta = mean(b); %mean of b
    e = X-G*b;      %error term
    rss = e'*e;     %residuals
    r=size(X,1)-rank(G); 
    sigma= sqrt(rss/r);
    
end
function [t] = glm_inference (c, b) %(c, beta, sigma, alpha)

    t=c*b/std(c*b);
end

function [X_mod, X_star,modelDataAllEffects] =glm_model(G,b,H,X, Gnoise)
    
    %%%Data model with all effects
    modelDataAllEffects= G(:,[1,2,3])*b([1,2,3],1);% + e;
    
    X_mod= G(:,[1,2,3,9,10])*b([1,2,3,9,10],1);% + e 
    X_star= X-(H*(b(4:11,1))) - Gnoise;              
end