function [varargout]=anlyModel(tspan,m,pullval)
% =================================================================== %
% function model to fit and predict GNSS time series
% Author: Sambuddha Dhar

% Input: tspan=time point; 
%          m=model parameter [ c1,c2,c3, tau1, tau2, lambda]
%         pullval='none' for normal; 1 for normalization
% Output: [dcal= modeled displacement
%           VR = displacements due to VR
%           AS = displacements due to AS
%           parscale= parameter scale if needed
% Example: [dc]=anlymodel(ts,m,'none');
%           dcal=c1+a*VR(Ar,t)+b*AS(A,t);

% See also: https://github.com/smbdhar/gpucb
%           email: sambuddha.dhar.p8@dc.tohoku.ac.jp

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
[a,b,c1,tom,K,A]=deal(m(1),m(2),m(3),m(4),m(5),m(6));
tref=max(tspan);
[strain,parstrain]=viscousflow(tom,K,tspan,tref,pullval);
[slip,parslip]=fltcreep(A,tspan,tref,pullval);
VR=a.*parstrain.*strain;  AS=b.*parslip.*slip; dcal=c1+VR+AS;
varargout{1}=dcal;  %;
varargout{2}=VR; %+(2*a)-c1
varargout{3}=AS; %+b-c1
varargout{4}=[a.*parstrain,b.*parslip,c1,tom,K,A];
end

function [strain,parstrain]=viscousflow(tom,K,tspan,tref,val)
n=3; m=3;
nn1=-1./(n-1); nn2=-1./(m-1);
% nn1=-(n); nn2=-(m);
% strain=(1+0.2)-((1+(tspan./tom)).^nn1) -(0.2.*(1+(tspan./(K))).^nn2);
strain=2-((1+(tspan./tom)).^nn1) -(((1+(tspan./K)).^nn2));
valm=interp1(tspan,strain,tref);
if strcmp(val,'none')
    parstrain=1;
else
    parstrain=val./valm;
end
end

function [slip,parslip]=fltcreep(A,tspan,tref,val)
% 0.0008 curve, 0.8 straight
G=30e3;
Lflt=90e3;
V0=0.08;%0.0599; % 
sn=300;
A=A.*sn;
slip=1-( 2.*A.*acoth(coth(1./(2.*A)).* exp((tspan.*2.*V0.*G)./(A.*Lflt))) );
% slip=0-( 2.*A.*acoth(coth(1./(2.*A)).* exp((tspan.*2.*V0.*G)./(A.*Lflt))) );
% slip=slip; %-(0.08.*tspan) (Lflt./G).*
% slip=slip-(tspan.*2.*V0.*G)./(A.*Lflt);
valm=interp1(tspan,slip,tref);
if strcmp(val,'none')
    parslip=1;
else
    parslip=val./valm;
end
end