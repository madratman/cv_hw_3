function [u, v] = LucasKanadeWithTemplateCorrection(It, It1, rect, I0, rect0)

thresh = 2;

[dun, dvn] = LucasKanade(It, It1, rect);
rect_n = rect + [dun, dvn, dun, dvn];

[u_corr, v_corr] = LucasKanadeWithInitialValue(I0, It1, rect0, (rect_n(1 : 2) - rect0(1 : 2))' );

% Update the outputs
u_corr = u_corr - (rect(1) - rect0(1));
v_corr = v_corr - (rect(2) - rect0(2));

if norm([u_corr, v_corr] - [dun, dvn]) <= thresh
    u = u_corr;
    v = v_corr;
else
    u = dun;
    v = dvn;
end
