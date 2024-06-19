% close all;
clear all;
c=clock; fprintf('RUN STARTED >>>>  Date: %d/%d/%d Time: %d:%d:%2.2f.\n',c(2),c(3),c(1),c(4),c(5),c(6));

if 0
    raw_data=dlmread('Well_Log_Example.prn');
    data(:,1)=raw_data(:,1);
    data(:,2)=raw_data(:,3);
    no_of_points=20;35;15;
else
    raw_data = load('G:\matlab\data\direct\gt\D2_011\4\tbc\acc_file.txt');
    raw_data1 = load('G:\matlab\data\direct\gt\D2_011\4\tbc\acc_file_knot.txt');
    
    acc = load('G:\matlab\data\direct\gt\D2_011\4\tbc\acc_file.txt');
    gyro = load('G:\matlab\data\direct\gt\D2_011\4\tbc\gyro_file.txt');
    figure,plot(acc(:,1), acc(:,2:4));hold on;plot(acc(:,1), acc(:,5:7));
    figure,plot(gyro(:,1), gyro(:,2:4));hold on;plot(gyro(:,1), gyro(:,5:7));
    
    delta = 40;
    figure;plot(raw_data(:,1), raw_data(:,2));hold on;plot(raw_data1(:,1), raw_data1(:,2),'or');
    plot(raw_data(:,1), raw_data(:,3) + delta);hold on;plot(raw_data1(:,1), raw_data1(:,3) + delta,'or');
    plot(raw_data(:,1), raw_data(:,4) + 2 * delta);hold on;plot(raw_data1(:,1), raw_data1(:,4) + 2 * delta,'or');
    grid on;
    figure,pcshow(raw_data(:,2:4),[1 0 0]);hold on;pcshow(raw_data(:,5:7),[0 0 1]);
    
    return
    
    raw_data2 = load('G:\matlab\data\direct\gt\D2_011\4\tbc\gyro_file.txt');
    raw_data0 = raw_data;
    raw_data = raw_data(1:40:end,:);
    [~,dist] = NormalizeVector(raw_data(:,2:4));
    if 1
        data=raw_data(:,1:4);
        data = [raw_data(:,1) dist raw_data(:,2:4)];
        no_of_points=round(size(data,1) / 10);
    else
        [xx,yy] = titanium;
        data(:,1)=xx;
        data(:,2)=yy;
%         pick = [1 5 11 21 27 29 31 33 35 40 45 49];
%         data = data(pick,:);
        no_of_points = 8;
    end
end


original_data=data;
original_size_data=size(data,1);
power = 3;
for j=1:original_size_data-no_of_points
    area=[];
    interpolated_dist=[];
    obj_func=[];
    a=original_size_data-no_of_points;
    x=data(:,1);
    y=data(:,2);
    xyz=data(:,3:5);
    for i=2:size(x,1)-1;
        xv=[x(i-1) x(i) x(i+1)];
        if 0
            yv=[y(i-1) y(i) y(i+1)];
        else
            yv=[xyz(i-1,1) xyz(i,1) xyz(i+1,1)];% x1 x2 x3
            yv2=[xyz(i-1,2) xyz(i,2) xyz(i+1,2)];% y1 y2 y3
            yv3=[xyz(i-1,3) xyz(i,3) xyz(i+1,3)];% z1 z2 z3
        end
        zv=[0 0 0];
        if 0
            xyz1 = [data([i-1:i+1],2) data([i-1:i+1],3) data([i-1:i+1],4)];
        else
            xyz1 = [xv' yv' [1 1 1]'];
        end
        AB = xyz1(1,:) - xyz1(3,:);  
        AC = xyz1(2,:) - xyz1(3,:);  
        crossProduct = cross(AB, AC);  
        triangleArea = norm(crossProduct) / 2;
        xyzv = [xyz(i-1,1:3); xyz(i,1:3); xyz(i+1,1:3)];
        if 1
            if 0
                area(i-1,1)=polyarea(xv,yv);
            else
                if 1
                    tempx =  polyarea(xv,yv);
                    tempy =  polyarea(xv,yv2);
                    tempz =  polyarea(xv,yv3);
                    area(i-1,1) = tempx^power + tempy^power + tempz^power;
                else
                   % 计算四面体体积
                   pt1 = [yv(1) yv2(1) yv3(1)];
                   pt2 = [yv(2) yv2(2) yv3(2)];
                   pt3 = [yv(3) yv2(3) yv3(3)];
                   plane = cross(pt1, pt2);
                   plane = plane./norm(plane);
                   dist = abs(dot(plane, pt3));
                   len1 = pt1 - [0 0 0];
                   len2 = pt2 - [0 0 0];
                   area_temp = norm(cross(pt1,pt2))/2;
                   volume = dist * area_temp / 3;
                   area(i-1,1) = volume;
                end
            end
        else
            area(i-1,1)=triangleArea;
        end
    end;
    obj_func=area;
    [minValue, rowsToDelete1] = min(obj_func(:,1));
    history_obj_func(j,1:2)=[original_size_data-j minValue];
    data(rowsToDelete1+1, :) = [];  % Get rid of min.
end;

if 0
    figure;
    plot(original_data(:,2),original_data(:,1),'b-');
    hold on;
    plot(data(:,2),data(:,1),'r:','LineWidth',2);
    hold on;
    plot(data(:,2),data(:,1),'ro','LineWidth',2);
%     set(gca,'YDir','Reverse');
    ylabel('Depth (ft)','FontName','Times New Roman','FontSize',14.0);
    xlabel('VP (ft/sec)','FontName','Times New Roman','FontSize',14.0);
    grid on;
else
    figure;
    subplot(2,2,1),plot(original_data(:,1),original_data(:,3),'b-x');
    hold on;
    plot(data(:,1),data(:,3),'r:','LineWidth',2);
    hold on;
    plot(data(:,1),data(:,3),'ro','LineWidth',2);
    control_point = 5 * ones(size(data,1),1);
    plot(data(:,1), control_point,'xm');
%     set(gca,'YDir','Reverse');
    ylabel('Depth (ft)','FontName','Times New Roman','FontSize',14.0);
    xlabel('VP (ft/sec)','FontName','Times New Roman','FontSize',14.0);
    title('x');
    grid on;
    subplot(2,2,2),plot(original_data(:,1),original_data(:,4),'b-x');
    hold on;
    plot(data(:,1),data(:,4),'r:','LineWidth',2);
    hold on;
    plot(data(:,1),data(:,4),'ro','LineWidth',2);
    control_point = 5 * ones(size(data,1),1);
    plot(data(:,1), control_point,'xm');
%     set(gca,'YDir','Reverse');
    ylabel('Depth (ft)','FontName','Times New Roman','FontSize',14.0);
    xlabel('VP (ft/sec)','FontName','Times New Roman','FontSize',14.0);
    title('y');
    grid on;
    subplot(2,2,3),plot(original_data(:,1),original_data(:,5),'b-x');
    hold on;
    plot(data(:,1),data(:,5),'r:','LineWidth',2);
    hold on;
    plot(data(:,1),data(:,5),'ro','LineWidth',2);
    control_point = 5 * ones(size(data,1),1);
    plot(data(:,1), control_point,'xm');
%     set(gca,'YDir','Reverse');
    ylabel('Depth (ft)','FontName','Times New Roman','FontSize',14.0);
    xlabel('VP (ft/sec)','FontName','Times New Roman','FontSize',14.0);
    title('z');
    grid on;
    subplot(2,2,4),plot(original_data(:,1),original_data(:,2),'b-x');
    hold on;
    plot(data(:,1),data(:,2),'r:','LineWidth',2);
    hold on;
    plot(data(:,1),data(:,2),'ro','LineWidth',2);
    control_point = 5 * ones(size(data,1),1);
    plot(data(:,1), control_point,'xm');
%     set(gca,'YDir','Reverse');
    ylabel('Depth (ft)','FontName','Times New Roman','FontSize',14.0);
    xlabel('VP (ft/sec)','FontName','Times New Roman','FontSize',14.0);
    title('dist');
    grid on;
    figure,plot3(raw_data0(:,2),raw_data0(:,3),raw_data0(:,4),'-');hold on;pcshow(raw_data0(:,2:4),[1 0 0]);hold on;pcshow(data(:,3:5),[0 0 1], 'MarkerSize', 50);
    figure,pcshow(raw_data0(:,2:4),[1 0 0]);hold on;pcshow(raw_data0(:,5:7),[0 0 1]);
end

figure;
plot(history_obj_func(:,1), history_obj_func(:,2),'b-','LineWidth',2);
ylabel('Error','FontName','Times New Roman','FontSize',14.0);
xlabel('Number of control points/knots','FontName','Times New Roman','FontSize',14.0);
grid on;