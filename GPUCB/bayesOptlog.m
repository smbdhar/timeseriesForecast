function [bestX,bestscore,Xr,y]=bayesOptlog(varargin)
% =================================================================== %
% Multivariate bayesian optimization GP-UCB method
% Author: Sambuddha Dhar
%
% Example; [bestX,bestscore,X,y]=bayes_opt(f,bounds,N,d,n,iter,fig)
% Input: f= loss function; log space
%           bounds = parameters bounds [m1min m2max; m2min m2max;...]
%           N= Initial sampling (nodal points)
%           n= reconnaissance points number (not related to sampling)
%           d= Dimension
%           iter=number of iteration;
%           fig= 'plot_true' for plotting the results
% Output: bestX=bestfit model parameters
%           bestscore= the value at best-fit parameters
%           Xr= the sampling point
%           Y= the value respect to sampling point

% See details: https://github.com/smbdhar/timeseriesForecast
%
% Note: Change wm1 wm2 in the plot section to see different parameters
%
% UCB method follows Srinivas et al. (2012 IEEE):https://doi.org/10.1109/TIT.2011.2182033

% Copyright (c) 2023 SAMBUDDHA DHAR
% 
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files 
% (the "Software"), to deal in the Software without restriction, including 
% without limitation the rights to use, copy, modify, merge, publish, 
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the 
% following conditions:
% 
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.
% =================================================================== %
f=varargin{1}; bounds=varargin{2}; N=varargin{3}; d=varargin{4}; n=varargin{5}; iter=varargin{6};
if nargin==7
    fig=varargin{7};
else
    fig='plot_false';
end
randunf=@(a,b,N,d) a+(b-a).*rand([N,d]);
% boundmat=[bounds, repmat([-1 1],size(bounds,1),1)];
% convert2Xi=@(Xmin,Xmax,Rmin,Rmax,Ri) (Xmin + (Ri-Rmin).*((Xmax-Xmin)./(Rmax-Rmin)) )';
% bounds=repmat([-1 1],size(bounds,1),1);
bounds=log10(bounds);

% initiate X and fX
X=randunf(bounds(:,1)',bounds(:,2)',N,d);
y=zeros(size(X,1),1);
for i=1:size(X,1)
%     XRi=convert2Xi(boundmat(:,1),boundmat(:,2),boundmat(:,3),boundmat(:,4),X(i,:)');
    y(i)=f(10.^X(i,:)); % <<Forward model>>
end

i0=[]; fxb0=[];
% if strcmp (fig, 'plot_true');f1=figure(91);f2=figure(92);end
for i=1:iter
    %initiation of test points
    Xtest=randunf(bounds(:,1)',bounds(:,2)',n,d);
    % calculaet mu and sigma of test points
    [mu_Xtest,sigma_Xtest]= return_musigma(X,y,Xtest);
    fxb=max(y);
    EI=UCB(mu_Xtest,sigma_Xtest,fxb,i);
    [~,bestEIid]=max(EI);
    pickid=randperm(length(bestEIid),1); bestEIid=bestEIid(pickid);
    bestXtest=Xtest(bestEIid,:);
%     bestXtestRi=convert2Xi(boundmat(:,1),boundmat(:,2),boundmat(:,3),boundmat(:,4),bestXtest');
    ybestXtest=f(10.^bestXtest); % <<Forward model>>
    % Update the search points
    X=[X;bestXtest];
    y=[y;ybestXtest];
    
    %------ plot ---------------------
    wm1=1; wm2=2;
    if strcmp (fig, 'plot_true') %fig==1
        figure(91);cla;hold on; box on;
%         Xr=convert2Xi(boundmat(:,1),boundmat(:,2),boundmat(:,3),boundmat(:,4),X');
        plot(X(:,wm1),X(:,wm2),'bo');%plot(10.^X(:,1),10.^X(:,2),'bo');
%         voronoi(10.^X(:,1),10.^X(:,2));
        bestpt=(unique(X(y==max(y),:),'rows'));
        plot(bestpt(:,wm1),bestpt(:,wm2),'r^');%plot(10.^bestpt(:,1),10.^bestpt(:,2),'r^');
        
%         xlim([10.^bounds(1,1) 10.^bounds(1,2)]); ylim([10.^bounds(2,1) 10.^bounds(2,2)])
        xlim([bounds(wm1,1) bounds(wm1,2)]); ylim([bounds(wm2,1) bounds(wm2,2)])
        title(sprintf('Iteration=%d;  fxb=%f',i,fxb))
        grid on; %view([15 45])
        
        fxb0=[fxb0;fxb]; i0=[i0;i];%ffxb(i)=fxb;
        figure(92);cla;%hold on;box on;
        plot(i0,fxb0,'r')
        axis([0 iter -0 1])        
        xlabel('iter');ylabel('error');
    end
    %------------------------------------
end
Xr=10.^X;
bestscore=unique(y(y==max(y),:));
bestX=10.^(unique(X(y==max(y),:),'rows'));
% bestX=convert2Xi(boundmat(:,1),boundmat(:,2),boundmat(:,3),boundmat(:,4),bestXr');
if strcmp (fig, 'plot_true')%fig==1
    figure(91);plot(bestX(:,wm1),bestX(:,wm2),'r^','MarkerSize',10);
end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
function [mu_Xtest,sigma_Xtest]= return_musigma(X,fX,Xtest)
% X=node points [sample x dimension]
% fX=f(X) where f() is loss function
% Xtest= test/quary points [sample x dimension]
nv=5e-14;    % noise variance; (nv)0.00000001
N=size(X,1); % Sample Number
fX=fX+ nv*rand(N,1);
K = kernelgen(X, X);
[L,flag] = chol(K + nv*eye(N)); L=L';
% Compute the mean at our test points.
Lk=(L\kernelgen(X, Xtest));
mu_Xtest = (Lk'*((L\fX)) )';%column
mu_Xtest =mu_Xtest';
% Compute the variance at our test points.
K_ = kernelgen(Xtest, Xtest);
s2 = diag(K_)' - sum(Lk.^2,1);
sigma_Xtest = sqrt(s2);
sigma_Xtest=sigma_Xtest';
end

function K=kernelgen(a,b)
% a,b = column vector; [sample x dimensions]
kp=10; % kernel parameter (K=20;kp=50)
K= exp((-0.5./kp).*( (sum(a.^2,2)+ sum(b.^2,2)'-2.*(a*b'))) );
% K=2.*K; %20
% n=1;
% K=( 1./sqrt(kp.*((2.*pi).^n)) ).*K;
end

function EI=UCB(mux,sigmax,fxb,iter)
% Calculate EI of test points
t=floor(iter./1)+0; D=2; del=0.1;
bt=sqrt(2.*log(D.*t.*t.*pi.*pi./(6.*del)))./5;
EI=mux +(bt.*sigmax); %UCB
% EI=mux +(0.5.*sigmax); %UCB
end
