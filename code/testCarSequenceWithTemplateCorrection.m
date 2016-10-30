load('../data/carseq.mat');
rect_init = [60, 117, 146, 152];
rects_all = zeros(size(frames, 3), 4);
rect_default = rect_init;
rect_temp_corr = rect_init;

for i = 1:size(frames, 3)-1
    [u, v] = LucasKanadeInverseCompositional(frames(:,:,i), frames(:,:,i+1), rect_default);
    rect_default = rect_default+[u,v,u,v];
    [u, v] = template_correction(frames(:, :, i), frames(:, :, i+1), rect_temp_corr, frames(:, :, 1), rect_init);
    rect_temp_corr = rect_temp_corr+[u,v,u,v];
    imshow(frames(:, :, i));
    hold on
    rectangle('Position', [rect_default(1), rect_default(2), rect_default(3)-rect_default(1), ...
        rect_default(4) - rect_default(2)], 'EdgeColor', 'r', 'LineWidth', 2);
    rectangle('Position', [rect_temp_corr(1), rect_temp_corr(2), rect_temp_corr(3)-rect_temp_corr(1), ...
        rect_temp_corr(4) - rect_temp_corr(2)], 'EdgeColor', 'g', 'LineWidth', 2);
    pause(0.01)
    axis off;
    if i==2|i==100|i==200|i==300|i==400
        saveas(gcf, strcat('car_templatecorr',int2str(i),'.png'));
    end
    rects_all(i, :) = rect_temp_corr;
end

save('../results/carseqrects-wcr.mat','rects_all');
close;