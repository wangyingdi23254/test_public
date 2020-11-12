%读入压力数据和放热率数据
clear
load('p_1.mat');
p(:,2)=p(:,2)/10;

%指定存放图像目录
path='E:\Experiment\Liuhui\LiuH\photo\1229\9_20141229_220610_C001H001S0001';
imdir=dir([path,'\*.bmp']);

myObj=VideoWriter('Case_1.avi');%初始化一个avi文件
myObj.FrameRate = 8;
open(myObj);
ind=[31:50:2913, 2913:1:2953];
imnum=length(ind);
for i=1:imnum
    imname=imdir(ind(i)).name;
    imind=str2double(imname(end-9:end-4));
    t=imind/288;
    p_int=interp1(p(:,1),p(:,2),t);
    im=imread([path,'\', imname]);
    im_pos=[0,0,0.5,1];
    subplot('Position',im_pos);
    imshow(im)
    p_pos=[0.6,0.59,0.37,0.24];
    subplot('Position',p_pos);
    plot(p(:,1),p(:,2));
    xlim([0 15]);
    ylim([0 22]);
    xlabel('Time after spark ignition / ms')
    ylabel('Pressure / MPa')
    set(gca,'Fontsize',8)
    hold on
    plot(t,p_int,'ro');
    hold off
    pmag_pos=[0.6,0.24,0.37,0.24];
    subplot('Position',pmag_pos);
    plot(p(:,1),p(:,2));
    xlabel('Time after spark ignition / ms')
    ylabel('Pressure / MPa')
    xlim([10 10.5]);
    ylim([0 22]);
    set(gca,'Fontsize',8)
    hold on
    plot(t,p_int,'ro');
    hold off
    frame = getframe(gcf);
    frame.cdata=frame.cdata(105:532,1:826,:);
    writeVideo(myObj,frame);
end
% movie2avi(M,'ani','compression','none','fps',8)
close(myObj);