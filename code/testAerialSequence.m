load('../data/aerialseq.mat');
for i=1:size(frames,3)-1
    i
    curr_mask = SubtractDominantMotion(frames(:,:,i),frames(:,:,i+1));
    img_fused = imfuse(frames(:,:,i), curr_mask);
    imshow(img_fused);
    pause(0.01)
    if (i==30)|(i==60)|(i==90)|(i==120)
        saveas(gcf, strcat('aerial',int2str(i), '.png'));
    end
end
close