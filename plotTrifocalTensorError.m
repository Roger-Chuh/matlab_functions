function plotTrifocalTensorError()
data = load('G:\matlab\data\direct\gt\D2_011\4\tbc\trifocal_tensor_result.txt');
assert(size(data,1)/4 == data(end,4) + 1)

err1 = data(1:4:end,6);
err2 = data(2:4:end,6);
err3 = data(3:4:end,6);
err4 = data(4:4:end,6);

err1 = err1(err1 > 0);
err2 = err2(err2 > 0);
err3 = err3(err3 > 0);
err4 = err4(err4 > 0);


figure,subplot(4,1,1);plot(err1);hold on;grid on;title(sprintf('camera [0 1 2], mean: %f, median: %f', mean(err1), median(err1)));
subplot(4,1,2);plot(err2);hold on;grid on;title(sprintf('camera [0 1 3], mean: %f, median: %f', mean(err2), median(err2)));
subplot(4,1,3);plot(err3);hold on;grid on;title(sprintf('camera [0 2 3], mean: %f, median: %f', mean(err3), median(err3)));
subplot(4,1,4);plot(err4);hold on;grid on;title(sprintf('camera [1 2 3], mean: %f, median: %f', mean(err4), median(err4)));

end