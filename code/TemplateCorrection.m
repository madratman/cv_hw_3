function [u, v] = template_correction(It, It1, rect, first_frame, rect_0)

thresh = 2;
[diff_u, diff_v] = LucasKanadeInverseCompositional(It, It1, rect);
rect_n = rect + [diff_u, diff_v, diff_u, diff_v];
uanvinthisbeautifulworld = (rect_n(1 : 2) - rect_0(1 : 2))'; 
p_init = [zeros(2,2) uanvinthisbeautifulworld];
[u_star, v_star] = LucasKanadeInverseCompositional_with_p_init(first_frame, It1, rect_0, p_init);
u_star = u_star - (rect(1)-rect_0(1));
v_star = v_star - (rect(2)-rect_0(2));
if norm([u_star, v_star]-[diff_u, diff_v])<=thresh
    [u v] = [u_star v_star];
else
    [u v] = [diff_u diff_v];
end