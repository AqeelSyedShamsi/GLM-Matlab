myfuns = stat;

i = 3000; %number of repetations
j = 1; %1000; %number of voxels
%c=ones(j,1);
alpha = 0.05;
loopItr = 10000;
betaMat=zeros (j,loopItr); %beta matrix
T=zeros (j,loopItr);       %T test matrix
%design matrix
for xx = 1:loopItr

    %creating gaussian signal and Design matrix
    
    [X,G] = myfuns.glm_simulate_data(i,j);
    

    %Parameters 
    [b, c, beta, e, rss, r, sigma] = myfuns.glm_estimate(X,G);
    %G = myfuns.glm_simulate_data(i,j);
    
    t = myfuns.glm_inference (c, sigma, b, G,i,rss);
    T (:, xx)= t;
    

end

p = 1-tcdf(T,i-1);
H1 = p < alpha;
disp(sum(H1)/10000)

D=6;
sigmaMat = eye(D);
detSigmaMat=det(sigmaMat);
invSigmaMat=(sigmaMat^(-1));

inputdata =rand(1,6);
x= random('normal',0,1,size(inputdata)); %RANDOM data


h = (1/((2*pi)^(D/2)))*(1/(detSigmaMat^(1/2)))*exp((-1/2)*((x-0)')*((x-0)*invSigmaMat));%Bishop 1.52