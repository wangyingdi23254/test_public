%%
path=uigetdir;	%选择图像文件目录
imfiles=dir([path '\*.bmp']);	%列出图像文件
imnum=length(imfiles)-1;	%图像数
imbg=imread([path '\bg.bmp']);	%读如背景图像
imsize=size(imbg);%	图像大小
datalength=imsize(1)*imsize(2);	
data=zeros(datalength,imnum,imsize(3)+1);
fc=8E3;	%截止频率
fs=3.6E5;	%拍摄频率
NFFT=2^(nextpow2(imnum)+0);	%FFT点数
f=fs/2*linspace(0,1,NFFT/2+1);	%FFT频率
f=(f/1000)';
yamp=zeros(datalength,NFFT/2+1,imsize(3)+1);	%振幅


%%
%%4层矩阵中，1-4层依次对应RGB分量和灰度级
for i=1:imnum
    imname=imfiles(i).name;
    im=imread([path '\' imname]);
    im=im-imbg;	%减背景
    for j=1:imsize(3)
        imComp=im(:,:,j);	%抽取颜色分量
        data(:,i,j)=reshape(imComp,datalength,1);	%将图像矩阵重整为向量
    end
    imgray=rgb2gray(im);	%转化为灰度图像
    data(:,i,end)=reshape(imgray,datalength,1);
end


%%
data_f=data;
imtime=(1:1:imnum)/fs;
for k=1:datalength;
    for m=1:imsize(3)+1
        pixs=data(k,:,m);
        pixs_f=avlfilter(imtime',pixs',fc,fs,'high');%高通滤波
        data_f(k,:,m)=pixs_f;
        y=fft(pixs_f,NFFT)/imnum;%FFT
        yamp(k,:,m)=2*abs(y(1:NFFT/2+1));
    end
end


%%
sumyamp=zeros(NFFT/2+1,imsize(3)+1);
for m=1:imsize(3)+1
    sumyamp(:,m)=sum(yamp(:,:,m))';
    makedirpath=[path, '\im_processed_',num2str(m)];
    mkdir(makedirpath);
end
figure,plot(f,sumyamp(:,1),'r',f,sumyamp(:,2),'g',f,sumyamp(:,3),'b',f,sumyamp(:,4),'k')
indf=1:NFFT/2+1;
figure,plot(indf,sumyamp(:,1),'r',indf,sumyamp(:,2),'g',indf,sumyamp(:,3),'b',indf,sumyamp(:,4),'k') 

%%
% ind10=21;
imf=zeros(imsize(1),imsize(2),imsize(3)+1);
imf2=imf;
makedirpath_0=[path, '\im_processed_',num2str(0)];
mkdir(makedirpath_0);
for n=1:NFFT/2+1
    for m=1:imsize(3)+1
        makedirpath=[path, '\im_processed_',num2str(m)];
%         imf0=reshape(yamp(:,ind10,m),imsize(1),imsize(2));
        imf(:,:,m)=reshape(yamp(:,n,m),imsize(1),imsize(2));    %将FFT幅值重整转化为矩阵
        imf2(:,:,m)=imf(:,:,m)*(255/max(max(imf(:,:,m))));
        imwrite(uint8(imf2(:,:,m)),[makedirpath, '\', 'im_processed_', num2str(n,'%03d'),'.bmp'],'bmp')
    end
    %     imf=imf*(255/max(max(imf0))); %以(1,0)模式为基增强显示
    imf(:,:,1:3)=imf(:,:,1:3)*(255/max(max(max(imf)))); %以各图像自身为基准增强显示
    imwrite(uint8(imf(:,:,1:3)),[makedirpath_0, '\', 'im_processed_', num2str(n,'%03d'),'.bmp'],'bmp')
end

%%
mkdir([path, '\im_hpf']);
for i=1:imnum
    imf=reshape(data_f(:,i,1:3),imsize(1),imsize(2),3); %将各像素高通滤波结果重新转化为矩阵
    imf=imf*5;	%增强显示
    imwrite(uint8(imf),[path, '\im_hpf\', '\im_hpf_', num2str(i,'%03d'),'.bmp'],'bmp')
end

%%
% P0=0;
% for j=1:datalength;
%    [S,F,T,P]=spectrogram(data(j,:),128,120,256,5E4);
%    P0=P+P0;
% end
%
% surf(T,F,10*log10(P),'edgecolor','none'); axis tight;
% view(0,90);
% xlabel('Time (Seconds)'); ylabel('Hz');



