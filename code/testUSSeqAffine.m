clear;
load('../data/usseq.mat');
load('../results/usseqrects.mat');
for i=1:size(frames,3)-1
    i;
    It = frames(:,:,i);
    It1 = frames(:,:,i+1);
    curr_rect = rects(i,:);
    next_rect = rects(i+1,:);
    curr_mask = SubtractDominantMotion(It(curr_rect(2):curr_rect(4),curr_rect(1):curr_rect(3)),...
                                       It1(next_rect(2):next_rect(4),next_rect(1):next_rect(3)));
    left_pad = zeros(size(It,1), round(curr_rect(1)));
    right_pad = zeros(size(It,1), round(size(It,2)-curr_rect(3)));
    up_pad = zeros(round(curr_rect(2)), size(curr_mask,2));
    down_pad = zeros(round(size(It,1)-size(up_pad,1)-size(curr_mask,1)), size(curr_mask,2));
    center_block = [up_pad;curr_mask;down_pad];

%     disp(['left pad ', num2str(size(left_pad))])
%     disp(['up_pad ', num2str(size(up_pad))])
%     disp(['curr_mask ', num2str(size(curr_mask))])
%     disp(['down_pad ', num2str(size(down_pad))])
%     disp(['right_pad ', num2str(size(right_pad))]
%     disp(['center_mask ', num2str(size(center_mask))])

    padded_mask = [left_pad center_block right_pad];

    img_fused = imfuse(It, padded_mask);
%     img_fused = imfuse(It, curr_mask);

    imshow(img_fused);
    pause(0.01)
    if (i==5)|(i==25)|(i==50)|(i==75)|(i==99)
        saveas(gcf, strcat('usseqaffine',int2str(i),'.png'));
    end
    axis off;
end
close