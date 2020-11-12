function p_filter=avlfilter(ca_intp,p_intp,fc,fs,method)
%ca_intp：x
%p_intp：y
%fc：截止频率
%fs：采样频率
p_filter=zeros(length(ca_intp),1);
n=ceil(2*fs/fc)+1;
p_add_L=ones(n,1)*p_intp(1);
p_add_R=ones(n,1)*p_intp(end);
p1=[p_add_L; p_intp; p_add_R];
p2=p1;
c=zeros(2*n+1,1);
h=zeros(2*n+1,1);
switch method
    case 'high'
        for j=1:2*n+1
            if j==n+1
                c(j)=1-2*fc/fs;
            else
                c(j)=-sin(2*pi*fc/fs*(j-(n+1)))/(pi*(j-(n+1)));
            end
            h(j)=0.5*(1+cos(2*pi*(j-(n+1))/(2*n+1)));
            k=h.*c;
        end
        for i=n+1:length(p1)-n
            p2(i)=sum(k.*(p1(i-n:i+n)));
        end
        p_filter=p2(n+1:length(p1)-n);
    case 'low'
        for j=1:2*n+1
            if j==n+1
                c(j)=2*fc/fs;
            else
                c(j)=sin(2*pi*fc/fs*(j-(n+1)))/(pi*(j-(n+1)));
            end
            h(j)=0.5*(1+cos(2*pi*(j-(n+1))/(2*n+1)));
            k=h.*c;
        end
        for i=n+1:length(p1)-n
            p2(i)=sum(k.*(p1(i-n:i+n)));
        end
        p_filter=p2(n+1:length(p1)-n);
end