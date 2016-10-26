function [u,v] = LucasKanadeBasis(It, It1, rect, bases)

% input - image at time t, image at t+1, rectangle (top left, bot right
% coordinates), bases 
% output - movement vector, [u,v] in the x- and y-directions.

It = double(It);
It1 = double(It1);

rect= round(rect);

[X,Y] = meshgrid(rect(1):rect(3), rect(2):rect(4));
sx = size(X,2);
sy = size(Y,1);
X = reshape(X,[sy*sx,1]);
Y = reshape(Y,[sy*sx,1]);

T = interp2(im2double(It),X,Y);   
[X,Y] = meshgrid(rect(1):rect(3), rect(2):rect(4));

T = reshape(T,[sy,sx]);       
patch = T;

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

p_arr_prev = zeros(6,1);
wrap_prev = [1+p_arr_prev(1) p_arr_prev(3) p_arr_prev(5); p_arr_prev(2) 1+p_arr_prev(4) p_arr_prev(6)];
steep_des_img = delta_matrix*jacobian;
steep_des_img_bases = reshape(steep_des_img, [size(patch),6]);

% eqn 46 from ref2 
for i=1:6
    for k=1:10
        steep_des_img_bases(:,:,i) = steep_des_img_bases(:,:,i) ...
                        - sum(sum(bases(:,:,k).*steep_des_img_bases(:,:,i),1),2)*bases(:,:,k);
    end
end

hess_all = zeros(6,6,sy,sx);
for i=1:sy
    for j=1:sx
      steep_des_img_bases_curr = reshape(steep_des_img_bases(i,j,:),[6,1]);              
      hess_all(:,:,i,j) = steep_des_img_bases_curr*steep_des_img_bases_curr';
    end
end        
hessian=sum(sum(hess_all,3),4);
hess_inv = pinv(hessian);
steep_des_img_bases = reshape(steep_des_img_bases, [size(patch,1)*size(patch,2), 6]);

while 1
    warped_coord = wrap_prev*[X(:)'; Y(:)'; ones([1,size(X(:))])];
    warped_patch = interp2(It1, reshape(warped_coord(1,:), size(X,1), size(X,2)),...
                           reshape(warped_coord(2,:), size(Y,1), size(Y,2)));
    error_image = warped_patch - patch;
    pre_final_matrix = hess_inv*steep_des_img_bases'*error_image(:);
    
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

    p_arr_prev = update_arr
    wrap_prev = [1+p_arr_prev(1) p_arr_prev(3) p_arr_prev(5); p_arr_prev(2) 1+p_arr_prev(4) p_arr_prev(6)];
    norm_to_check=norm(delta_p_matrix(5:6));
    if norm(norm_to_check)<0.001
        break;
    end
end    

u = round(p_arr_prev(5),5);
v = round(p_arr_prev(6),5);