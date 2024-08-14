function test_exposure()
inputDir = 'G:\matlab\data\direct\gt\D2_011\4\tbc\exposure11_dark\SLAMRecord';
dirInfo = dir(fullfile(inputDir, '*'));

info = {};
ids = zeros(length(dirInfo)-2,1);
for i = 3 : length(dirInfo)
name = dirInfo(i).name;
idx = find(name == '_');
info{i-2,1}.id = str2double(name(1:idx(1)-1));
ids(i-2,1) = info{i-2,1}.id;
info{i-2,1}.exposure = str2double(name(idx(1)+1:idx(2)-1));
info{i-2,1}.gain = str2double(name(idx(2)+1:idx(3)-1));
end
[~,idx] = sort(ids);
info2  = cell(length(info),1);


for i = 1 : length(ids)
   info2{i,1} = info{idx(i),1};    
end

camera_index = [1 2 4 5];
exposure_gain_mat = zeros(length(info2),4);
for i = 3 : length(dirInfo)
   path = fullfile(inputDir, dirInfo(idx(i-2)+2).name); 
   ind = find(dirInfo(idx(i-2)+2).name == '_');
   check(i-2,1) = str2double(dirInfo(idx(i-2)+2).name(1:ind(1)-1));
   pathInfo = dir(fullfile(path, 'Camera*'));
   camid  = 1;
   imgs = cell(4,1);
   for folder_id = camera_index
      json = loadjson(fullfile(path, pathInfo(folder_id).name,'data.json')); 
      if folder_id == 1
          info2{i-2}.exposure_json = json.Sequence.Frameset.Frame(1).exposure_time;
          info2{i-2}.gain_json = json.Sequence.Frameset.Frame(1).gain;
          exposure_gain_mat(i-2,:) = [info2{i-2}.exposure info2{i-2}.gain info2{i-2}.exposure_json info2{i-2}.gain_json];
      end
      imgInfo = dir(fullfile(path, pathInfo(folder_id).name,'images','*.bmp'));
      imgs{camid,1} = imread(fullfile(imgInfo(end).folder, imgInfo(end).name));
      camid = camid + 1;
   end
   img_big = [imgs{2,1} imgs{3,1};imgs{1,1} imgs{4,1}];
   imwrite(img_big, fullfile('G:\matlab\data\direct\gt\D2_011\4\tbc\exposure11_dark', sprintf('%03d_%07d_%04d_%07d_%04d.png',i-2, exposure_gain_mat(i-2,:))));
end

end