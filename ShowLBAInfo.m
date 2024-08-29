function ShowLBAInfo()

data = load('G:\matlab\data\direct\gt\D2_011\4\tbc\lba_res.txt');

all_timestamps = unique(data(:,1));

err_info = [];
for i = 1 : length(all_timestamps)
   id1 = find(data(:,1) == all_timestamps(i)); 
%    id2 = find(data(:,2) == all_timestamps(i)); 
   temp = data(id1,:);
   weighted_mean = sum(temp(:,3).*temp(:,5)) ./ sum(temp(:,5));
   raw_mean = sum(temp(:,4).*temp(:,5)) ./ sum(temp(:,5));
   
   weighted_mean_inlier = sum(temp(:,6).*temp(:,8)) ./ sum(temp(:,8));
   raw_mean_inlier = sum(temp(:,7).*temp(:,8)) ./ sum(temp(:,8));
   err_info = [err_info; [all_timestamps(i) weighted_mean raw_mean weighted_mean_inlier raw_mean_inlier sum(temp(:,8))/sum(temp(:,5)) sum(temp(end,5)) sum(temp(end,8))]];
end


figure,plot(repmat(err_info(:,1), 1, 7), err_info(:,2:8));grid on;legend('mean of weighted err', 'mean of raw err', 'mean of weighted err (inliers only)', 'mean of raw err (inliers only)','inlier ratio','all vm in cur fid','inlier vm in cur fid');title('inlier threshold: 0.5 pixel');
end