Turnover_T=[];
Turnover_tao=[];
for i=1:22
    x0=1000./T(:,i);
    y0=(tid_LLNL(:,i));
    h= abs(diff([x0(2), x0(1)]));
    y1=gradient(y0, h);
    [xi,yi,iout,jout]=intersections(x0,y1,[1000/1000, 1000/600],[0,0]);
    yi=interp1(x0,y0,xi);
    if ~isempty(xi)
        Turnover_T=[Turnover_T,xi];
        Turnover_tao=[Turnover_tao,yi];
    end
end
semilogy(1000./T(:,1),tid_LLNL,'color',[0.7,0.7,0.7] );
hold on
% plot(Turnover_T(1,:),Turnover_tao(1,:),'color',[246 139 31]/255,'linewidth',1)
% plot(Turnover_T(2,:),Turnover_tao(2,:),'color',[246 139 31]/255,'linewidth',1)

[~,ind_p_bottom]=size(Turnover_T);
ntc_top_T=Turnover_T(1,1):1E-3:Turnover_T(2,1);
ntc_top_tao=interp1(1000./T(:,1),tid_LLNL(:,1),ntc_top_T);
ntc_bottom_T=Turnover_T(1,end):1E-3:Turnover_T(2,end);
ntc_bottom_tao=interp1(1000./T(:,ind_p_bottom),tid_LLNL(:,ind_p_bottom),ntc_bottom_T);

fillx=[Turnover_T(1,:),ntc_bottom_T,flip(Turnover_T(2,:)),flip(ntc_top_T)];
filly=[Turnover_tao(1,:),ntc_bottom_tao,flip(Turnover_tao(2,:)),flip(ntc_top_tao)];
fill(fillx,filly,[246 139 31]/255,'facealpha',0.1,'edgecolor',[246 139 31]/255)
hold off