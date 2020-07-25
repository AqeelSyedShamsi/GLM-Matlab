function func = moreVoxelBOLDstat
   
    func.glm_simulate_data=@glm_simulate_data;
    func.glm_estimate=@glm_estimate;
    func.glm_inference=@glm_inference;
    func.glm_model=@glm_model;

end
function [X, G, H, Gnoise, resp, hrf] = glm_simulate_data(i,j,TR,HR,BR, gamma,res,GMean,GStd)

    G=ones(i,1);%G = eye(i,j); %design matrix
    Gnoise = GMean + GStd.*randn(i,1);%Gaussian Noise with mean 1 and standard deviation 2.

    hrf = spm_hrf(1/res);
    resp=[];
    hrf1=hrf(1:30);
    for respi = 1:i/60
        resp=vertcat(resp, hrf1);
    end
    
    
    %n1=1:1:i;
    %n=linspace (0,0.001,300) ;
    n=0.01:0.01:3;
    t=n*TR;
    hrFreq=sin(2*pi*HR*t);
    brFreq=sin(2*pi*BR*t);

    
    inputdata =rand(i,6);
    H= (random('normal',0,GStd,size(inputdata))); %RANDOM data
    H =[H hrFreq'];
    H =[H brFreq'];
    
   
    %resp= resp + Gnoise + (H*gamma');
    %X= resp+ Gnoise + (H*gamma');
    var1=vertcat(resp,zeros(length(resp),1));
    var2=vertcat(zeros(length(resp),1),resp);
    X=[(var1+Gnoise+(H*gamma')) (var2+Gnoise+(H*gamma'))];
    fftCondA= fft(X(:,1));
    fftCondB= fft(X(:,2));
    
    absFftCondA = abs(fftCondA);
    absFftCondB = abs(fftCondB);

    magFftCondA = max(absFftCondA);
    magFftCondB = max(absFftCondB);
    rho=magFftCondA / magFftCondB;

    
    %X=X+Gnoise;
    G= [G X H];
    
    %X=[X ()];
    
    
    
end
function [b, beta, e, rss, r, sigma, inveG] =glm_estimate(X,G)
    inveG=inv(G'*G);
    %b = pinv(G)*X;
    
    b= (inveG*G')*X;
   
    beta = mean(b);
    e = X-G*b;
    rss = e'*e;
    r=size(X,1)-rank(G);
    sigma= sqrt(rss/r);
    
end
function [t] = glm_inference (c, b) %(c, beta, sigma, alpha)

    t=c*b/std(c*b) ;
    
end
