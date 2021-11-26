function [MinRMSE,NoClasses]=Umethod(x,y,fig)
%
%
%  Description: Adapted from this function applies the L-method (Salvador and Chan, 2005) for the estimation of the
%              appropriate number of classes on an evaluation graph.
%
%  Inputs:
%    x: x-coordinates of the evaluation graph (number of classes)
%    y: y-coordinates of the evaluation graph (value of a distance/similarity/error/quality)
%    fig: 'on' or 'off'
%  Outputs:
%    rmse :  The total root mean squared error
%
%    NoClasses: the appropriate number of clusters based on the intersection
%               of two best-fit lines to the left and right side around the
%               'knee', a critical location on the evaluation graph denoted
%               with a green square marker.
%
%
%  Example:
%           close all, clc, clear all
%           x=[1:16]';
%           y=[44.51; 33.71; 26.11; 21.05; 17.19; 14.32; 12.50;	...
%              10.51; 7.63; 4.47; 2.47; 2.02; 2.00; 2.00; 2.00; 2.00];
%           [rmse,NoClasses]=Lmethod(x,y)
%
%
%  Figures:
%       upper left: the evaluation graph
%       upper right: all the possible best-fit lines (green & red pairs)
%       lower left: the minimization of RMSE (a green square marker denotes the 'knee')
%       lower right: the best-fit lines correnspoding to the minimum RMSE
%
%
%  *** Further details about this function can be found in the paper :
%       Salvador, S., Chan, P., 2005. Learning States and Rules for Detecting
%       Anomalies in Time Series. Applied Intelligence 23(3), 241-255.
%
%  Version 2.00 (23/10/2012) <----
%  Version 1.00 (17/02/2012)


% ************************************************************************************************************
% ************************************************************************************************************

% If you use this code for academic use, it would be appreciated to cite the following publication(s) in you work:
%
% *** [1] A. Zagouras, R.H. Inman, C.F.M. Coimbra, On the determination of coherent solar microclimates for utility planning and operations,
% Solar Energy, Volume 102, April 2014, Pages 173-188,
% ISSN 0038-092X, http://dx.doi.org/10.1016/j.solener.2014.01.021.
%
% *** [2] A. Zagouras, A. Kazantzidis, E. Nikitidou, A.A. Argiriou, Determination of measuring sites for solar irradiance, based on cluster analysis of satellite-derived cloud estimations,
% Solar Energy, Volume 97, November 2013, Pages 1-11,
% ISSN 0038-092X, http://dx.doi.org/10.1016/j.solener.2013.08.005.
%
%
% Any questions or comments are mostly welcome; please contact me at my email address below.
%
%
%  Copyright (c) 2012, Athanassios Zagouras, University of California, San Diego (UCSD)
%  All rights reserved.
%  Email: azagouras@eng.ucsd.edu; thzagouras@gmail.com;
% ************************************************************************************************************
% ************************************************************************************************************

if nargin<3
    fig='off'
end


ffont=10;
axfsize=10;
MS=12;
switch fig
    case 'on'
        subplot(2,2,1);
        plot(x,y,'b.-','MarkerSize', MS); xlim([0 max(x)+2]);
        title('Evaluation Graph','FontSize',ffont,'FontWeight','bold')
        xlabel('Number of classes','FontSize',ffont);
        ylabel('Evaluation Metric','FontSize',ffont);
        
        
        hold on;
end

t=0;
d=length(x);
for b=2:d-4
    for c=b+2:d-2; t=t+1;
        %-------------------------------&
        p1=1:b;
        p2=b+1:c;
        p3=c+1:d;
        
        X1=x(p1);
        Y1=y(p1);
        
        X2=x(p2);
        Y2=y(p2);
        
        X3 = x(p3);
        Y3 = y(p3);
        
        coef_1st=polyfit(X1,Y1,1);
        yfit1st = coef_1st(1)*X1 + coef_1st(2);
        rmse1st=sqrt(sum((y(p1)-yfit1st).^2)/numel(y(p1)));
        
        coef_2nd=polyfit(X2,Y2,1);
        yfit2nd = coef_2nd(1)*X2 + coef_2nd(2);
        rmse2nd=sqrt(sum((y(p2)-yfit2nd).^2)/numel(y(p2)));
        
        coef_3rd=polyfit(X3,Y3,1);
        yfit3rd = coef_3rd(1)*X3 + coef_3rd(2);
        rmse3rd=sqrt(sum((y(p3)-yfit3rd).^2)/numel(y(p3)));
        
        %         plot(x(p1),y(p1),'r-')
        %         plot(x(p2),y(p2),'g-')
        %         plot(x(p3),y(p3),'b-')
        
        b_factor=x(b);
        c_factor=x(c);
        d_factor=x(end);
        
        rmse_tot(t)=(b_factor-1)./(d_factor-1)*rmse1st + ...
            (c_factor-b_factor)./(d_factor-1)*rmse2nd +...
            (d_factor-c_factor)./(d_factor-1)*rmse3rd;
        
        B(t)=b;
        C(t)=c;
    end 
end

[MinRMSE,pos]=min(rmse_tot); BestPos=pos+1; 
NoClasses=[B(BestPos) C(BestPos)];

switch fig
    case 'on'
        subplot(2,2,3); plot3(B,C,rmse_tot,'k.', 'MarkerSize',3); hold on ;
        plot3(B(BestPos),C(BestPos),rmse_tot(BestPos),'r.:', 'MarkerSize',MS);
        title('RMSE','FontSize',ffont,'FontWeight','bold')
        xlabel('Number of points in first line','FontSize',ffont);
        ylabel('Number of points in 2nd line','FontSize',ffont);
        zlabel('RMSE_c','FontSize',ffont);
        set(gca,'fontsize',axfsize);
        
        subplot(2,2,4);
        scatter(B,C,10,rmse_tot); hold on
        colorbar
        plot(B(BestPos),C(BestPos),'r.:', 'MarkerSize',MS);
        xlabel('Number of points in first line','FontSize',ffont);
        ylabel('Number of points in 2nd line','FontSize',ffont);
        set(gca,'fontsize',axfsize);
        title('RMSE','FontSize',ffont,'FontWeight','bold')
end

switch fig
    case 'on'
p1=1:B(BestPos);
p2=B(BestPos)+1:C(BestPos);
p3=C(BestPos)+1:d;

X1=x(p1);
Y1=y(p1);

X2=x(p2);
Y2=y(p2);

X3 = x(p3);
Y3 = y(p3);

coef_1st=polyfit(X1,Y1,1);
yfit1st = coef_1st(1)*X1 + coef_1st(2);
rmse1st=sqrt(sum((y(p1)-yfit1st).^2)/numel(y(p1)));

coef_2nd=polyfit(X2,Y2,1);
yfit2nd = coef_2nd(1)*X2 + coef_2nd(2);
rmse2nd=sqrt(sum((y(p2)-yfit2nd).^2)/numel(y(p2)));

coef_3rd=polyfit(X3,Y3,1);
yfit3rd = coef_3rd(1)*X3 + coef_3rd(2);
rmse3rd=sqrt(sum((y(p3)-yfit3rd).^2)/numel(y(p3)));

        subplot(2,2,2)
        plot(x(p1),y(p1),'ro'); hold on;
        plot(x(p2),y(p2),'go')
        plot(x(p3),y(p3),'bo')
        plot(x(p1),yfit1st,'k-','linewidth',2); hold on;
        plot(x(p2),yfit2nd,'k-','linewidth',2)
        plot(x(p3),yfit3rd,'k-','linewidth',2)
        
        title('Best Fit Lines','FontSize',ffont,'FontWeight','bold')
        xlabel('Number of classes','FontSize',ffont);
        ylabel('Evaluation Metric','FontSize',ffont);
end
end
