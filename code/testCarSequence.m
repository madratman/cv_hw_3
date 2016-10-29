load('../data/carseq.mat');
rect=[60, 117, 146, 152];
rects = zeros(size(frames,3),4);
rects(1,:)=rect;

for i=1:size(frames,3)-1
    It=frames(:,:,i);
    It1=frames(:,:,i+1);
    It=double(It);
    It1=double(It1);
    [u,v] = LucasKanadeInverseCompositional(It, It1, rect);
    rect = rect+[u v u v];
    rects(i+1,:)=rect;
    imagesc(It)
    hold on; 
    rectangle('Position',[rect(1), rect(2), abs(rect(3)-rect(1)), abs(rect(2)-rect(4))],...
        'EdgeColor','r','LineWidth',3);
    pause(0.01);
    if i==2|i==100|i==200|i==300|i==400
        saveas(gcf, strcat('car',int2str(i), '.png'));
    end
    axis off;
end
save('../results/carseqrect.mat','rects');
close;