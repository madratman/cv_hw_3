function M = LucasKanadeAffine(It, It1)

% input - image at time t, image at t+1 
% output - M affine transformation matrix

It = double(It);
It1 = double(It1);

patch = It;
[X,Y] = meshgrid(1:size(It,2), 1:size(It,1));
[delta_x, delta_y] = gradient(patch);
steep_des_img = zeros(numel(patch),6);

for i=1:numel(patch)
    steep_des_img(i,:) = [delta_x(i) delta_y(i)]*[X(i) 0 Y(i) 0 1 0; 0 X(i) 0 Y(i) 0 1];
end

p_arr_prev = zeros(6,1);
warp_prev = [1+p_arr_prev(1) p_arr_prev(3) p_arr_prev(5); p_arr_prev(2) 1+p_arr_prev(4) p_arr_prev(6)];

while 1
    warped_coord = warp_prev*[X(:)'; Y(:)'; ones([1,size(X(:))])];
    warped_patch = interp2(It1, reshape(warped_coord(1,:), size(X,1), size(X,2)),...
                                reshape(warped_coord(2,:), size(Y,1), size(Y,2)));
    error_image = warped_patch - patch;
    error_image = error_image(:);
    error_image(isnan(error_image)) = 0;
    steep_des_img(isnan(error_image),:) = [];
    hess_curr = steep_des_img'*steep_des_img;
    hess_inv = pinv(hess_curr);
    
    pre_final_matrix = hess_inv*steep_des_img'*error_image;
    delta_p_matrix = [-pre_final_matrix(1)-pre_final_matrix(1)*pre_final_matrix(4)+pre_final_matrix(2)*pre_final_matrix(3);
                      -pre_final_matrix(2);
                      -pre_final_matrix(3);
                      -pre_final_matrix(4)-pre_final_matrix(1)*pre_final_matrix(4)+pre_final_matrix(2)*pre_final_matrix(3);
                      -pre_final_matrix(5)-pre_final_matrix(4)*pre_final_matrix(5)+pre_final_matrix(3)*pre_final_matrix(6);
                      -pre_final_matrix(6)-pre_final_matrix(1)*pre_final_matrix(6)+pre_final_matrix(2)*pre_final_matrix(5)];
                  
    den_factor = (1+pre_final_matrix(1))*(1+pre_final_matrix(4)) - pre_final_matrix(2)*pre_final_matrix(3);
    delta_p_matrix = delta_p_matrix/den_factor;
    
    update_arr = [p_arr_prev(1)+delta_p_matrix(1)+p_arr_prev(1)*delta_p_matrix(1)+p_arr_prev(3)*delta_p_matrix(2);
                  p_arr_prev(2)+delta_p_matrix(2)+p_arr_prev(2)*delta_p_matrix(1)+p_arr_prev(4)*delta_p_matrix(2);
                  p_arr_prev(3)+delta_p_matrix(3)+p_arr_prev(1)*delta_p_matrix(3)+p_arr_prev(3)*delta_p_matrix(4);
                  p_arr_prev(4)+delta_p_matrix(4)+p_arr_prev(2)*delta_p_matrix(3)+p_arr_prev(4)*delta_p_matrix(4);
                  p_arr_prev(5)+delta_p_matrix(5)+p_arr_prev(1)*delta_p_matrix(5)+p_arr_prev(3)*delta_p_matrix(6);
                  p_arr_prev(6)+delta_p_matrix(6)+p_arr_prev(2)*delta_p_matrix(5)+p_arr_prev(4)*delta_p_matrix(6)];

    p_arr_prev = update_arr;
    warp_prev = [1+p_arr_prev(1) p_arr_prev(3) p_arr_prev(5); p_arr_prev(2) 1+p_arr_prev(4) p_arr_prev(6)];

    if norm(pre_final_matrix)<0.1
        break;
    end
end    

M = [warp_prev; 0 0 1];