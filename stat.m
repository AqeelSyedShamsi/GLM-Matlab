function func = stat
   
    func.glm_simulate_data=@glm_simulate_data;
    func.glm_estimate=@glm_estimate;
    func.glm_inference=@glm_inference;
end


function [X, G] = glm_simulate_data(i,j)
    inputdata =rand(i,j);
    X= random('normal',0,1,size(inputdata)); %RANDOM data
    G=ones(i,1);%G = eye(i,j); %design matrix
    
end

function [b, c, beta, e, rss, r, sigma] =glm_estimate(X,G)
    b = pinv(G)*X;
    d=-1*(b<0);
    dd=1*(b>0);
    c = d+dd;
    beta = mean(b);
    e = X-G*b;
    rss = e'*e;
    r=size(X,1)-rank(G);
    sigma= sqrt(rss/r);
    %inveG=inv(G'*G);
end
function [t,Cb] = glm_inference (c, sigma, b, G,i,rss) %(c, beta, sigma, alpha)

    Cb = (sigma^2)*(inv(G'*G));
    %t=(c'*b)/(sigma*sqrt(c*inveG*c'));
    %t = (c'*b)./(sqrt((rss/(i-rank(G))))*sqrt(c*inveG*c'));
    t = (1*b)./(sqrt((rss/(i-rank(G))))*sqrt(inv(G'*G)));
    
    
end

