load('../data/sylvseq.mat');
load('../data/sylvbases.mat');
rect_bases = [102, 62, 156, 108];
rect_inv_comp = [102, 62, 156, 108];
rects_all = zeros(size(frames,3),4);
rects_all(1,:)=rect_bases;
for i=1:length(frames)-1
     i
    [u_bases,v_bases]=LucasKanadeBasis(frames(:,:,i),frames(:,:,i+1),rect_bases,bases);
    [u_inv_comp,v_inv_comp]=LucasKanadeInverseCompositional(frames(:,:,i),frames(:,:,i+1),rect_inv_comp);
    rect_bases = rect_bases + [ u_bases v_bases u_bases v_bases];
    rect_inv_comp = rect_inv_comp + [ u_inv_comp v_inv_comp u_inv_comp v_inv_comp];
    rects_all(i+1,:)=rect_bases;
    imshow(frames(:,:,i))
    hold on
    rectangle('Position',[rect_inv_comp(1),rect_inv_comp(2),...
        rect_inv_comp(3)-rect_inv_comp(1),rect_inv_comp(4)-rect_inv_comp(2)],'EdgeColor','b'); 
    rectangle('Position',[rect_bases(1),rect_bases(2),...
        rect_bases(3)-rect_bases(1),rect_bases(4)-rect_bases(2)],'EdgeColor','r'); 
    pause(0.01);     
    if (i==2)|(i==200)|(i==300)|(i==350)|(i==400)
        saveas(gcf, strcat('sylv',int2str(i), '.png'));
    end
 end  
save('../results/slyvseqrect.mat','rects');