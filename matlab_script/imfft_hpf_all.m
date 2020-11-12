%%
path=uigetdir;	%ѡ��ͼ���ļ�Ŀ¼
imfiles=dir([path '\*.bmp']);	%�г�ͼ���ļ�
imnum=length(imfiles)-1;	%ͼ����
imbg=imread([path '\bg.bmp']);	%���米��ͼ��
imsize=size(imbg);%	ͼ���С
datalength=imsize(1)*imsize(2);	
data=zeros(datalength,imnum,imsize(3)+1);
fc=8E3;	%��ֹƵ��
fs=3.6E5;	%����Ƶ��
NFFT=2^(nextpow2(imnum)+0);	%FFT����
f=fs/2*linspace(0,1,NFFT/2+1);	%FFTƵ��
f=(f/1000)';
yamp=zeros(datalength,NFFT/2+1,imsize(3)+1);	%���


%%
%%4������У�1-4�����ζ�ӦRGB�����ͻҶȼ�
for i=1:imnum
    imname=imfiles(i).name;
    im=imread([path '\' imname]);
    im=im-imbg;	%������
    for j=1:imsize(3)
        imComp=im(:,:,j);	%��ȡ��ɫ����
        data(:,i,j)=reshape(imComp,datalength,1);	%��ͼ���������Ϊ����
    end
    imgray=rgb2gray(im);	%ת��Ϊ�Ҷ�ͼ��
    data(:,i,end)=reshape(imgray,datalength,1);
end


%%
data_f=data;
imtime=(1:1:imnum)/fs;
for k=1:datalength;
    for m=1:imsize(3)+1
        pixs=data(k,:,m);
        pixs_f=avlfilter(imtime',pixs',fc,fs,'high');%��ͨ�˲�
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
        imf(:,:,m)=reshape(yamp(:,n,m),imsize(1),imsize(2));    %��FFT��ֵ����ת��Ϊ����
        imf2(:,:,m)=imf(:,:,m)*(255/max(max(imf(:,:,m))));
        imwrite(uint8(imf2(:,:,m)),[makedirpath, '\', 'im_processed_', num2str(n,'%03d'),'.bmp'],'bmp')
    end
    %     imf=imf*(255/max(max(imf0))); %��(1,0)ģʽΪ����ǿ��ʾ
    imf(:,:,1:3)=imf(:,:,1:3)*(255/max(max(max(imf)))); %�Ը�ͼ������Ϊ��׼��ǿ��ʾ
    imwrite(uint8(imf(:,:,1:3)),[makedirpath_0, '\', 'im_processed_', num2str(n,'%03d'),'.bmp'],'bmp')
end

%%
mkdir([path, '\im_hpf']);
for i=1:imnum
    imf=reshape(data_f(:,i,1:3),imsize(1),imsize(2),3); %�������ظ�ͨ�˲��������ת��Ϊ����
    imf=imf*5;	%��ǿ��ʾ
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



