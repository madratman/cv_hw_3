function [u, v] = LucasKanadeInverseCompositional_with_p_init(It, It1, rect, p)

It = double(It);
It1 = double(It1);

patch = It(rect(2):rect(4), rect(1):rect(3));
[X,Y] = meshgrid(rect(1):rect(3), rect(2):rect(4));
[delta_x, delta_y] = gradient(patch);

jacobian = zeros(2*numel(patch), 6);

for i=1:numel(patch)
    jacobian(2*i-1,:) = [0 0 0 0 1 0];
    jacobian(2*i, :) = [0 0 0 0 0 1];
end

delta_matrix = zeros(numel(patch), 2*numel(patch));

for i=1:numel(patch)
    delta_matrix(i, 2*i-1) = delta_x(i);
    delta_matrix(i, 2*i) = delta_y(i);
end

% initialize with given p
p_arr_prev = p; 
wrap_prev = [1+p_arr_prev(1) p_arr_prev(3) p_arr_prev(5); p_arr_prev(2) 1+p_arr_prev(4) p_arr_prev(6)];
steep_des_img = delta_matrix*jacobian;
hessian = steep_des_img'*steep_des_img;
hess_inv = pinv(hessian);

while 1
    warped_coord = wrap_prev*[X(:)'; Y(:)'; ones([1,size(X(:))])];
    warped_patch = interp2(It1, reshape(warped_coord(1,:), size(X,1), size(X,2)),...
                           reshape(warped_coord(2,:), size(Y,1), size(Y,2)));
    error_image = warped_patch - patch;
    pre_final_matrix = hess_inv*steep_des_img'*error_image(:);
    
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
    wrap_prev = [1+p_arr_prev(1) p_arr_prev(3) p_arr_prev(5); p_arr_prev(2) 1+p_arr_prev(4) p_arr_prev(6)];

    if norm(delta_p_matrix)<0.1
        break;
    end
end    

u = p_arr_prev(5);
v = p_arr_prev(6);