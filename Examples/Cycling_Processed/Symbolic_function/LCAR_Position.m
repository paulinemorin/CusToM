function out1 = LCAR_Position(in1,in2,in3)
%LCAR_POSITION
%    OUT1 = LCAR_POSITION(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    25-Jan-2021 14:40:27

R5cut1_1 = in3(37);
R5cut1_2 = in3(40);
R5cut1_3 = in3(43);
R5cut2_1 = in3(38);
R5cut2_2 = in3(41);
R5cut2_3 = in3(44);
R5cut3_1 = in3(39);
R5cut3_2 = in3(42);
R5cut3_3 = in3(45);
p5cut1 = in2(13);
p5cut2 = in2(14);
p5cut3 = in2(15);
q38 = in1(38,:);
q39 = in1(39,:);
q40 = in1(40,:);
q41 = in1(41,:);
q42 = in1(42,:);
t2 = cos(q38);
t3 = cos(q39);
t4 = cos(q40);
t5 = cos(q41);
t6 = cos(q42);
t7 = sin(q38);
t8 = sin(q39);
t9 = sin(q40);
t10 = sin(q41);
t11 = sin(q42);
t12 = R5cut1_1.*t2;
t13 = R5cut1_2.*t3;
t14 = R5cut1_3.*t2;
t15 = R5cut2_1.*t2;
t16 = R5cut2_2.*t3;
t17 = R5cut2_3.*t2;
t18 = R5cut3_1.*t2;
t19 = R5cut3_2.*t3;
t20 = R5cut3_3.*t2;
t21 = R5cut1_1.*t7;
t22 = R5cut1_2.*t8;
t23 = R5cut1_3.*t7;
t24 = R5cut2_1.*t7;
t25 = R5cut2_2.*t8;
t26 = R5cut2_3.*t7;
t27 = R5cut3_1.*t7;
t28 = R5cut3_2.*t8;
t29 = R5cut3_3.*t7;
t30 = -t23;
t31 = -t26;
t32 = -t29;
t33 = t14+t21;
t34 = t17+t24;
t35 = t20+t27;
t36 = t12+t30;
t37 = t15+t31;
t38 = t18+t32;
t39 = t9.*t33;
t40 = t9.*t34;
t41 = t9.*t35;
t42 = t3.*t36;
t43 = t3.*t37;
t44 = t3.*t38;
t45 = t8.*t36;
t46 = t8.*t37;
t47 = t8.*t38;
t48 = -t45;
t49 = -t46;
t50 = -t47;
t51 = t22+t42;
t52 = t25+t43;
t53 = t28+t44;
t54 = t13+t48;
t55 = t16+t49;
t56 = t19+t50;
t57 = t4.*t51;
t58 = t4.*t52;
t59 = t4.*t53;
t60 = -t57;
t61 = -t58;
t62 = -t59;
t63 = t39+t60;
t64 = t40+t61;
t65 = t41+t62;
out1 = [R5cut1_2.*(-3.333342798044145e-1)+p5cut1-t13.*2.800843105023983e-1+t45.*2.800843105023983e-1-t10.*t54.*1.847715543158421e-2+t5.*t63.*1.847715543158421e-2-t11.*(t4.*t33+t9.*t51).*9.035329006044679e-2-t6.*(t5.*t54+t10.*t63).*9.035329006044679e-2;R5cut2_2.*(-3.333342798044145e-1)+p5cut2-t16.*2.800843105023983e-1+t46.*2.800843105023983e-1-t10.*t55.*1.847715543158421e-2+t5.*t64.*1.847715543158421e-2-t11.*(t4.*t34+t9.*t52).*9.035329006044679e-2-t6.*(t5.*t55+t10.*t64).*9.035329006044679e-2;R5cut3_2.*(-3.333342798044145e-1)+p5cut3-t19.*2.800843105023983e-1+t47.*2.800843105023983e-1-t10.*t56.*1.847715543158421e-2+t5.*t65.*1.847715543158421e-2-t11.*(t4.*t35+t9.*t53).*9.035329006044679e-2-t6.*(t5.*t56+t10.*t65).*9.035329006044679e-2];
