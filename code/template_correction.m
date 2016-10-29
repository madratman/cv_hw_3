function [u, v] = template_correction(frame_n, frame_nplusone, rect_given, I0, rect0)

epsil = 2
[delta_u_n, delta_v_n] = LucasKanadeInverseCompositional(frame_n, frame_nplusone, rect_given);
rect_n = rect_given + [delta_u_n, delta_v_n, delta_u_n, delta_v_n];
u_and_v = (rect_n(1 : 2) - rect0(1 : 2))';
p_init = [zeros(2,2) u_and_v];
[u_star, v_star] = LucasKanadeInverseCompositional_with_p_init(I0, frame_nplusone, rect0, p_init);

u_star = u_star - (rect_given(1) - rect0(1));
v_star = v_star - (rect_given(2) - rect0(2));

if norm([u_star, v_star] - [delta_u_n, delta_v_n]) <= epsil
    u = u_star;
    v = v_star;
else
    u = delta_u_n;
    v = delta_v_n;
end