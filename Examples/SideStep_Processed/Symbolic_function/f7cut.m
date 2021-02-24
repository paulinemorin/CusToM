function [R7cut,p7cut] = f7cut(in1,in2,in3)
%F7CUT
%    [R7CUT,P7CUT] = F7CUT(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    17-Feb-2021 14:52:16

R6cut1_1 = in3(46);
R6cut1_2 = in3(49);
R6cut1_3 = in3(52);
R6cut2_1 = in3(47);
R6cut2_2 = in3(50);
R6cut2_3 = in3(53);
R6cut3_1 = in3(48);
R6cut3_2 = in3(51);
R6cut3_3 = in3(54);
p6cut1 = in2(16);
p6cut2 = in2(17);
p6cut3 = in2(18);
q16 = in1(13,:);
q17 = in1(14,:);
q43 = in1(37,:);
q44 = in1(38,:);
q45 = in1(39,:);
q46 = in1(40,:);
t2 = cos(q16);
t3 = cos(q17);
t4 = cos(q43);
t5 = cos(q44);
t6 = cos(q45);
t7 = cos(q46);
t8 = sin(q16);
t9 = sin(q17);
t10 = sin(q43);
t11 = sin(q44);
t12 = sin(q45);
t13 = sin(q46);
t14 = R6cut1_1.*t2;
t15 = R6cut1_2.*t2;
t16 = R6cut1_3.*t3;
t17 = R6cut2_1.*t2;
t18 = R6cut2_2.*t2;
t19 = R6cut2_3.*t3;
t20 = R6cut3_1.*t2;
t21 = R6cut3_2.*t2;
t22 = R6cut3_3.*t3;
t23 = R6cut1_1.*t8;
t24 = R6cut1_2.*t8;
t25 = R6cut1_3.*t9;
t26 = R6cut2_1.*t8;
t27 = R6cut2_2.*t8;
t28 = R6cut2_3.*t9;
t29 = R6cut3_1.*t8;
t30 = R6cut3_2.*t8;
t31 = R6cut3_3.*t9;
t32 = -t23;
t33 = -t26;
t34 = -t29;
t35 = t14+t24;
t36 = t17+t27;
t37 = t20+t30;
t38 = t15+t32;
t39 = t18+t33;
t40 = t21+t34;
t41 = t3.*t35;
t42 = t3.*t36;
t43 = t3.*t37;
t44 = t9.*t35;
t45 = t9.*t36;
t46 = t9.*t37;
t47 = -t41;
t48 = -t42;
t49 = t5.*t38;
t50 = -t43;
t51 = t5.*t39;
t52 = t5.*t40;
t53 = t11.*t38;
t54 = t11.*t39;
t55 = t11.*t40;
t59 = t16+t44;
t60 = t19+t45;
t61 = t22+t46;
t56 = -t53;
t57 = -t54;
t58 = -t55;
t62 = t25+t47;
t63 = t28+t48;
t64 = t31+t50;
t65 = t4.*t59;
t66 = t4.*t60;
t67 = t4.*t61;
t68 = t10.*t59;
t69 = t10.*t60;
t70 = t10.*t61;
t71 = t4.*t62;
t72 = t4.*t63;
t73 = t4.*t64;
t74 = t10.*t62;
t75 = t10.*t63;
t76 = t10.*t64;
t77 = -t74;
t78 = -t75;
t79 = -t76;
t80 = t68+t71;
t81 = t69+t72;
t82 = t70+t73;
t83 = t65+t77;
t84 = t66+t78;
t85 = t67+t79;
t86 = t6.*t80;
t87 = t6.*t81;
t88 = t6.*t82;
t89 = t5.*t83;
t90 = t5.*t84;
t91 = t5.*t85;
t92 = t11.*t83;
t93 = t11.*t84;
t94 = t11.*t85;
t95 = t49+t92;
t96 = t51+t93;
t97 = t52+t94;
t98 = t56+t89;
t99 = t57+t90;
t100 = t58+t91;
t101 = -t12.*(t53-t89);
t102 = -t12.*(t54-t90);
t103 = -t12.*(t55-t91);
t104 = t86+t101;
t105 = t87+t102;
t106 = t88+t103;
R7cut = reshape([t13.*t95-t7.*t104,t13.*t96-t7.*t105,t13.*t97-t7.*t106,t7.*t95+t13.*t104,t7.*t96+t13.*t105,t7.*t97+t13.*t106,-t12.*t80-t6.*(t53-t89),-t12.*t81-t6.*(t54-t90),-t12.*t82-t6.*(t55-t91)],[3,3]);
if nargout > 1
    p7cut = [p6cut1+t15.*4.816e-3+t16.*7.047784313725489e-2-t23.*4.816e-3-t25.*4.803577777777777e-2+t41.*4.803577777777777e-2+t44.*7.047784313725489e-2-t49.*3.1992e-1-t92.*3.1992e-1;p6cut2+t18.*4.816e-3+t19.*7.047784313725489e-2-t26.*4.816e-3-t28.*4.803577777777777e-2+t42.*4.803577777777777e-2+t45.*7.047784313725489e-2-t51.*3.1992e-1-t93.*3.1992e-1;p6cut3+t21.*4.816e-3+t22.*7.047784313725489e-2-t29.*4.816e-3-t31.*4.803577777777777e-2+t43.*4.803577777777777e-2+t46.*7.047784313725489e-2-t52.*3.1992e-1-t94.*3.1992e-1];
end
