function func = BOLDstat
   
    func.glm_simulate_data=@glm_simulate_data;
    func.glm_estimate=@glm_estimate;
    func.glm_inference=@glm_inference;
    func.glm_model=@glm_model;

end


function [X, G, H, Gnoise, resp, stim, hrf] = glm_simulate_data(i,j,TR,HR,BR, gamma,res,GMean,GStd)
    %inputdata =rand(i,j);
    %X= random('normal',0,1,size(inputdata)); %RANDOM data
    G=ones(i,1);%G = eye(i,j); %design matrix
    Gnoise = GMean + GStd.*randn(i,1);%Gaussian Noise with mean 1 and standard deviation 2.
    
%     hrf = spm_hrf(1/res);
%     % create a set of stimuli occurring at specific times % in the time series. First make an array of zeros, % then put ones where you want to have activity occur 
     stim=zeros(i,1);
%     %stim( 1*res:5*res:end ) =1;
     stim( 10:30:end ) =1;
%     resp = conv(stim, hrf);
%     resp = resp(1:length(stim));
%     resp = resp / max(resp);

    hrf = spm_hrf(1/res);
    resp=[];
    hrf1=hrf(1:30);
    for respi = 1:i/30
        resp=vertcat(resp, hrf1);
    end
    
    %n1=1:1:i;
    %n=linspace (0,0.001,300) ;
    n=0.01:0.01:3;
    t=n*TR;
    hrFreq=sin(2*pi*HR*t);
    brFreq=sin(2*pi*BR*t);

    
    inputdata =rand(i,6);
    H= (random('normal',0,1,size(inputdata))); %RANDOM data
    H =[H hrFreq'];
    H =[H brFreq'];
    
    G= [G resp H];
    
    X= resp+ Gnoise + (H*gamma');
    
end

function [b, beta, e, rss, r, sigma, inveG] =glm_estimate(X,G)
    b = pinv(G)*X;
   
    beta = mean(b);
    e = X-G*b;
    rss = e'*e;
    r=size(X,1)-rank(G);
    sigma= sqrt(rss/r);
    inveG=inv(G'*G);
end
function [t,Cb] = glm_inference (c, sigma, b, G,i,rss,inveG) %(c, beta, sigma, alpha)

    Cb = (sigma^2)*(inveG);
    %t=(c'*b)./(sigma*sqrt(c*inveG*c'));
    %t=(c'*b)/((Cb)^2);
    %t = (c'*b)./(sqrt((rss/(i-rank(G))))*sqrt(c*inveG*c'));
    t = (c*b)/(((rss/(i-rank(G)))*(c*inveG*c'))^0.5);
    %t = (1*b)./(sqrt((rss/(i-rank(G))))*sqrt(inv(G'*G)));
    
end
function [X_mod, X_star,modelDataAllEffects] =glm_model(G,e,b,H,gamma,X, Gnoise)
    
    %%%Data model with all effects
    modelDataAllEffects= G(:,[1,2])*b([1,2],1);% + e;
    
    %litH = H(:,[5,6,7,8]);
    %litGamma=gamma(:,[5,6,7,8]);
    X_mod= G(:,[1,2])*b([1,2],1);% + e
    X_star= X-(H*(b(3:10,1))) - Gnoise;
    %X_star = X-(H*((mean(H,1))')) - Gnoise ;
end