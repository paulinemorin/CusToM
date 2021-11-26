function [ind]=Best3Lines(x,y)
% ii=1
% x=dCOPcut{ii}(:,1);
% y=fct_moyenne_glissante(dCOPcut{ii}(:,2),25);
[~,ind] = min(abs(y-((max(y)-min(y))/2+min(y))));

c=ind;
p12=1:c;
X12 =x(p12);
Y12 =y(p12);
%plot(X12,Y12)

p3 = (floor(length(x)/2)+1):length(x);
X3 = x(p3);
Y3 = y(p3);

[~,No,~,yfit1st,~,yfit2nd]=Lmethod(X12,Y12,'on');
% b = NoClasses(1,1);
b = find(x==No);
rmse1st=sqrt(sum((y(1:b)-yfit1st).^2)/numel(y(1:b)));
rmse2nd=sqrt(sum((y(b+1:c)-yfit2nd).^2)/numel(y(b+1:c)));

plot(x(1:b),y(1:b))
plot(x(b+1:c),y(b+1:c))

coef_3rd=polyfit(X3,Y3,1);
yfit3rd = coef_3rd(1)*X3 + coef_3rd(2);
rmse3rd=sqrt(sum((y(p3)-yfit3rd).^2)/numel(y(p3)));

plot(X3,yfit3rd,'k-')
xlim('auto')

b_factor=x(b);
c_factor=x(c);
d_factor=x(end);

rmse_tot(1)=(b_factor-1)./(d_factor-1)*rmse1st + ...
    (c_factor-b_factor)./(d_factor-1)*rmse2nd +...
    (d_factor-c_factor)./(d_factor-1)*rmse3rd;

 B(1)=b;
 C(1)=c;

cmpt=1;
Vrmse(cmpt)=1;

while Vrmse>0.05
    cmpt=cmpt+1;
    if mod(cmpt,2)==0
        b=B(cmpt-1);
        p1=1:b;
        X1 =x(p1);
        Y1 =y(p1);
        
        p23 = (b+1):length(x);
        X23 = x(p23);
        Y23 = y(p23);
        
        [~,No,~,yfit2nd,~,yfit3rd]=Lmethod(X23,Y23,'off');
        % b = NoClasses(1,1);
        c = find(x==No);
        rmse2nd=sqrt(sum((y(b+1:c)-yfit2nd).^2)/numel(y(b:c)));
        rmse3rd=sqrt(sum((y(c+1:end)-yfit3rd).^2)/numel(y(c+1:end)));
        
        coef_1st=polyfit(X1,Y1,1);
        yfit1st = coef_1st(1)*X1 + coef_1st(2);
        rmse1st=sqrt(sum((y(p1)-yfit1st).^2)/numel(y(p1)));
        
        b_factor=x(b);
        c_factor=x(c);
        d_factor=x(end);

        rmse_tot(cmpt)=(b_factor-1)./(d_factor-1)*rmse1st + ...
            (c_factor-b_factor)./(d_factor-1)*rmse2nd +...
            (d_factor-c_factor)./(d_factor-1)*rmse3rd;
    else
        c=C(cmpt-1);
        p12=1:c;
        X12 =x(p12);
        Y12 =y(p12);
        
        p3 = (c+1):length(x);
        X3 = x(p3);
        Y3 = y(p3);
        
        [~,No,~,yfit1st,~,yfit2nd]=Lmethod(X12,Y12,'off');
        b = find(x==No);
        rmse1st=sqrt(sum((y(1:b)-yfit1st).^2)/numel(y(1:b)));
        rmse2nd=sqrt(sum((y(b+1:c)-yfit2nd).^2)/numel(y(b+1:c)));
        
        coef_3rd=polyfit(X3,Y3,1);
        yfit3rd = coef_3rd(1)*X3 + coef_3rd(2);
        rmse3rd=sqrt(sum((y(p3)-yfit3rd).^2)/numel(y(p3)));
        
        b_factor=x(b);
        c_factor=x(c);
        d_factor=x(end);
        
        rmse_tot(cmpt)=(b_factor-1)./(d_factor-1)*rmse1st + ...
            (c_factor-b_factor)./(d_factor-1)*rmse2nd +...
            (d_factor-c_factor)./(d_factor-1)*rmse3rd;
    end
    B(cmpt)=b;
    C(cmpt)=c;
    Vrmse(cmpt)=abs((rmse_tot(cmpt)-rmse_tot(cmpt-1))/rmse_tot(cmpt-1));
end
ind=[B(cmpt), C(cmpt)];
end