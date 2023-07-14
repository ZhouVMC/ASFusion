function J = SSIF(I,G,radius,Epsilon,kappa,scale)
padMethod = 'symmetric'; 
patchSize = 2*radius + 1;
h = ones(patchSize)/patchSize/patchSize;    

mu = imfilter(I, h, padMethod); % patch mean of I
nu = imfilter(G, h,padMethod); % patch mean of G
phi = imfilter(I.*G, h,padMethod) - mu.*nu; % patch covariance 

varSigma = max(0, imfilter(G.*G, h,padMethod) - nu.*nu); %patch var of G

a = phi./(varSigma + Epsilon);
Beta = (a + sign(phi).*sqrt(a.^2 + 4*kappa*Epsilon./(varSigma + Epsilon)))/2;

%weight calculation
w = varSigma./(scale*mean(varSigma(:)));
w = 1./(1+w.^2);
normalizeFactor = imfilter(w,h,padMethod); 

%final output
A = imfilter(Beta.*w, h,padMethod);
B = imfilter((mu - Beta.*nu).*w, h,padMethod);
J = (G.*A + B)./normalizeFactor;
end

