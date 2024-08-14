function WriteGif(inputDir, gifName,delay, inputDir2)
dirInfo = dir(fullfile(inputDir, '*.png'));
dirInfo2 = dir(fullfile(inputDir2, '*.png'));

if length(dirInfo) == 0
    dirInfo = dir(fullfile(inputDir, '*.bmp'));
end
if length(dirInfo2) == 0
    dirInfo2 = dir(fullfile(inputDir2, '*.bmp'));
end

nImages = length(dirInfo);
filename = fullfile(inputDir, gifName); % Specify the output file name
[a,b,c] = size(imread(fullfile(inputDir, dirInfo(end).name)));

mean_gray = zeros(nImages,1);
mean_gray_each = zeros(nImages,4);
mean_gray2 = zeros(nImages,1);
mean_gray_each2 = zeros(nImages,4);
for idx = 1 : nImages
    %     img = imresize(imread(fullfile(inputDir, dirInfo(idx).name)),[360,640]);
    %     img = imresize(imread(fullfile(inputDir, dirInfo(idx).name)),[922,1914]);
    img = imresize(imread(fullfile(inputDir, dirInfo(idx).name)),[1040,1920]);
    
    img = imresize(imread(fullfile(inputDir, dirInfo(idx).name)),[480,640]);
    img = imresize(imread(fullfile(inputDir, dirInfo(idx).name)),[480,1280]);
    img = imresize(imread(fullfile(inputDir, dirInfo(idx).name)),round([640,1706])./1);
    img = imresize(imread(fullfile(inputDir, dirInfo(idx).name)),round([960,1280])./1);
    img = imresize(imread(fullfile(inputDir, dirInfo(idx).name)),round([740,920])./1);
    img0 = imread(fullfile(inputDir, dirInfo(idx).name));
    img2 = double(img0(1:0.5*size(img0,1),1:0.5*size(img0,2)));
    img3 = double(img0(1:0.5*size(img0,1),0.5*size(img0,2)+1:end));
    img1 = double(img0(0.5*size(img0,1)+1:end,1:0.5*size(img0,2)));
    img4 = double(img0(0.5*size(img0,1)+1:end,0.5*size(img0,2)+1:end));
    
    mean_gray_each(idx,1) = sum(img1(:)) / size(img1,1) / size(img1,2);
    mean_gray_each(idx,2) = sum(img2(:)) / size(img2,1) / size(img2,2);
    mean_gray_each(idx,3) = sum(img3(:)) / size(img3,1) / size(img3,2);
    mean_gray_each(idx,4) = sum(img4(:)) / size(img4,1) / size(img4,2);
    
    mean_gray(idx,1) = sum(sum(double(img0))) / size(img0,1) / size(img0,2);
    img = imresize(img0, 0.5);
    img = insertText(img, [50, 50], sprintf('level: %d', str2double(dirInfo(idx).name(1:3))), 'FontSize', 20, 'TextColor', 'red');
    
    if length(dirInfo2) ~= 0
        img0 = imread(fullfile(inputDir2, dirInfo(idx).name));
        img2 = double(img0(1:0.5*size(img0,1),1:0.5*size(img0,2)));
        img3 = double(img0(1:0.5*size(img0,1),0.5*size(img0,2)+1:end));
        img1 = double(img0(0.5*size(img0,1)+1:end,1:0.5*size(img0,2)));
        img4 = double(img0(0.5*size(img0,1)+1:end,0.5*size(img0,2)+1:end));
        
        mean_gray_each2(idx,1) = sum(img1(:)) / size(img1,1) / size(img1,2);
        mean_gray_each2(idx,2) = sum(img2(:)) / size(img2,1) / size(img2,2);
        mean_gray_each2(idx,3) = sum(img3(:)) / size(img3,1) / size(img3,2);
        mean_gray_each2(idx,4) = sum(img4(:)) / size(img4,1) / size(img4,2);
        
        mean_gray2(idx,1) = sum(sum(double(img0))) / size(img0,1) / size(img0,2);
        img_ = imresize(img0, 0.5);
        img_ = insertText(img_, [50, 50], sprintf('level: %d', str2double(dirInfo(idx).name(1:3))), 'FontSize', 20, 'TextColor', 'red');
        img = [img img_];
    end
    
    
    %     img = img(:,1:1706/2,:);
    if(ndims(img) ~= 3)
        img = cat(3, img, img, img);
    end
    % %     img = imresize(img(21:386,143:583,:), 2);
    [A,map] = rgb2ind(img,256);
    if idx == 1
        for j = 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',delay);
        end
    else
        %         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',2);
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',delay);
    end
end
end