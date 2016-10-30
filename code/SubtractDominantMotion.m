function mask = SubtractDominantMotion(image1, image2)

% input - image1 and image2 form the input image pair
% output - mask is a binary image of the same size

image1 = double(image1);
image2 = double(image2);

M = LucasKanadeAffine(image1, image2);

%fill with nan
image_warped = warpH(image1, M, [size(image2, 1), size(image2, 2)], NaN);
image_warped(isnan(image_warped)) = 0;
mask = abs(image2 - image_warped);
mask = mask>37.5;
% thresh = graythresh(mask);
% mask = im2bw(mask, thresh);
SE = strel('disk', 8);
mask = imdilate(mask, SE);
mask = imerode(mask, SE);
mask = bwareaopen(mask, 50);
% imagesc(mask)
