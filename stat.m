function func = stat
   
    func.glm_simulate_data=@glm_simulate_data;
    func.glm_estimate=@glm_estimate;
    func.glm_inference=@glm_inference;
end


function [X, G, H] = glm_simulate_data(i,j,TR,HR,BR, gamma)
    inputdata =rand(i,j);
    X= random('normal',0,1,size(inputdata)); %RANDOM data
    G=ones(i,1);%G = eye(i,j); %design matrix
    
    
    
    
    %n1=1:1:i;
    %n=linspace (0,0.001,300) ;
    n=0.01:0.01:3;
    t=n*TR;
    hrFreq=sin(2*pi*HR*t);
    brFreq=sin(2*pi*BR*t);

    
    inputdata =rand(i,6);
    H= (random('normal',1.5,11,size(inputdata))); %RANDOM data
    H =[H hrFreq'];
    H =[H brFreq'];
    
    G= [G H];
    
    X= X + H*gamma';
    
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
   
    t = (c*b)/(((rss/(i-rank(G)))*(c*inveG*c'))^0.5);
    
end
