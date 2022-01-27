function [Rcut,pcut] = fcut(in1)
%FCUT
%    [RCUT,PCUT] = FCUT(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    18-Jan-2022 15:41:18

q1 = in1(1,:);
q2 = in1(2,:);
q3 = in1(3,:);
q4 = in1(4,:);
q5 = in1(5,:);
q6 = in1(6,:);
q7 = in1(7,:);
q8 = in1(8,:);
q9 = in1(9,:);
q10 = in1(10,:);
q11 = in1(11,:);
q12 = in1(12,:);
q13 = in1(13,:);
q17 = in1(17,:);
q18 = in1(18,:);
q19 = in1(19,:);
q20 = in1(20,:);
q21 = in1(21,:);
q24 = in1(24,:);
q25 = in1(25,:);
q26 = in1(26,:);
q27 = in1(27,:);
q28 = in1(28,:);
q31 = in1(31,:);
q32 = in1(32,:);
q38 = in1(38,:);
q39 = in1(39,:);
q46 = in1(45,:);
q47 = in1(46,:);
q48 = in1(47,:);
q49 = in1(48,:);
q50 = in1(49,:);
t2 = cos(q1);
t3 = cos(q2);
t4 = cos(q3);
t5 = cos(q4);
t6 = cos(q5);
t7 = cos(q6);
t8 = cos(q7);
t9 = cos(q8);
t10 = cos(q9);
t11 = cos(q10);
t12 = cos(q11);
t13 = cos(q12);
t14 = cos(q13);
t15 = cos(q17);
t16 = cos(q18);
t17 = cos(q19);
t18 = cos(q20);
t19 = cos(q21);
t20 = cos(q24);
t21 = cos(q25);
t22 = cos(q26);
t23 = cos(q27);
t24 = cos(q28);
t25 = cos(q31);
t26 = cos(q32);
t27 = cos(q38);
t28 = cos(q39);
t29 = cos(q49);
t30 = cos(q50);
t31 = sin(q1);
t32 = sin(q2);
t33 = sin(q3);
t34 = sin(q4);
t35 = sin(q5);
t36 = sin(q6);
t37 = sin(q7);
t38 = sin(q8);
t39 = sin(q9);
t40 = sin(q10);
t41 = sin(q11);
t42 = sin(q12);
t43 = sin(q13);
t44 = sin(q17);
t45 = sin(q18);
t46 = sin(q19);
t47 = sin(q20);
t48 = sin(q21);
t49 = sin(q24);
t50 = sin(q25);
t51 = sin(q26);
t52 = sin(q27);
t53 = sin(q28);
t54 = sin(q31);
t55 = sin(q32);
t56 = sin(q38);
t57 = sin(q39);
t58 = sin(q49);
t59 = sin(q50);
t60 = t2.*t29;
t61 = t2.*t58;
t62 = t29.*t31;
t63 = t4.*t59;
t64 = t16.*t59;
t65 = t21.*t59;
t66 = t31.*t58;
t67 = t33.*t59;
t68 = t45.*t59;
t69 = t50.*t59;
t71 = t2.*t3.*t30;
t72 = t2.*t15.*t30;
t73 = t2.*t20.*t30;
t74 = t4.*t29.*t30;
t75 = t16.*t29.*t30;
t76 = t21.*t29.*t30;
t77 = t2.*t30.*t32;
t78 = t3.*t30.*t31;
t79 = t2.*t30.*t44;
t80 = t15.*t30.*t31;
t81 = t2.*t30.*t49;
t82 = t20.*t30.*t31;
t84 = t4.*t30.*t58;
t85 = t29.*t30.*t33;
t86 = t16.*t30.*t58;
t87 = t29.*t30.*t45;
t88 = t21.*t30.*t58;
t89 = t29.*t30.*t50;
t90 = t30.*t31.*t32;
t91 = t30.*t31.*t44;
t92 = t30.*t31.*t49;
t95 = t30.*t33.*t58;
t96 = t30.*t45.*t58;
t97 = t30.*t50.*t58;
t195 = t59.*8.558498637340108e-2;
t204 = t30.*t31.*4.81415548350381e-2;
t214 = t30.*t31.*5.349061648337568e-2;
t216 = t2.*t30.*8.558498637340109e-3;
t217 = t29.*t30.*8.558498637340108e-2;
t219 = t30.*t58.*8.558498637340108e-2;
t236 = t2.*t30.*5.349061648337565e-3;
t70 = t59.*t66;
t83 = t59.*t60;
t93 = t59.*t61;
t94 = t59.*t62;
t99 = -t74;
t100 = -t75;
t101 = -t76;
t103 = -t90;
t104 = -t91;
t105 = -t92;
t106 = -t95;
t107 = -t96;
t108 = -t97;
t111 = t77+t78;
t112 = t79+t80;
t113 = t81+t82;
t202 = t60.*4.81415548350381e-2;
t203 = t61.*4.81415548350381e-2;
t205 = t67.*1.277082454163948e-1;
t212 = t60.*5.349061648337568e-2;
t213 = t61.*5.349061648337568e-2;
t215 = -t204;
t218 = t62.*8.558498637340109e-3;
t220 = t66.*8.558498637340109e-3;
t227 = t66.*5.349061648337565e-3;
t228 = t85.*1.277082454163948e-1;
t229 = t95.*1.277082454163948e-1;
t234 = -t216;
t238 = t62.*5.349061648337565e-3;
t98 = -t70;
t102 = -t83;
t109 = t61+t94;
t110 = t62+t93;
t128 = t71+t103;
t129 = t72+t104;
t130 = t73+t105;
t131 = t4.*t111;
t132 = t16.*t112;
t133 = t21.*t113;
t134 = t33.*t111;
t135 = t45.*t112;
t136 = t50.*t113;
t230 = t94.*4.81415548350381e-2;
t231 = -t212;
t232 = t70.*4.81415548350381e-2;
t233 = -t213;
t235 = -t218;
t237 = -t220;
t239 = -t229;
t240 = t94.*5.349061648337568e-2;
t241 = t70.*5.349061648337568e-2;
t243 = t83.*8.558498637340109e-3;
t244 = t93.*8.558498637340109e-3;
t247 = t83.*5.349061648337565e-3;
t248 = t93.*5.349061648337565e-3;
t114 = t60+t98;
t115 = t66+t102;
t116 = t3.*t109;
t117 = t3.*t110;
t118 = t15.*t109;
t119 = t15.*t110;
t120 = t20.*t109;
t121 = t20.*t110;
t122 = t32.*t109;
t123 = t32.*t110;
t124 = t44.*t109;
t125 = t44.*t110;
t126 = t49.*t109;
t127 = t49.*t110;
t152 = -t131;
t153 = t5.*t128;
t154 = -t132;
t155 = t17.*t129;
t156 = -t133;
t157 = t22.*t130;
t158 = t34.*t128;
t165 = t63+t134;
t166 = t64+t135;
t167 = t65+t136;
t242 = -t232;
t245 = -t240;
t246 = -t244;
t249 = -t247;
t278 = t131.*1.277082454163948e-1;
t137 = t3.*t114;
t138 = t3.*t115;
t139 = t15.*t114;
t140 = t15.*t115;
t141 = t20.*t114;
t142 = t20.*t115;
t143 = t32.*t114;
t144 = t32.*t115;
t145 = -t123;
t146 = t44.*t114;
t147 = t44.*t115;
t148 = -t125;
t149 = t49.*t114;
t150 = t49.*t115;
t151 = -t127;
t162 = -t153;
t163 = -t155;
t164 = -t157;
t168 = t67+t152;
t169 = t68+t154;
t170 = t69+t156;
t171 = t5.*t165;
t172 = t34.*t165;
t173 = t46.*t166;
t174 = t51.*t167;
t279 = -t278;
t301 = t153.*1.686712675310875e-2;
t159 = -t144;
t160 = -t147;
t161 = -t150;
t175 = t6.*t168;
t176 = t18.*t169;
t177 = t23.*t170;
t178 = t35.*t168;
t179 = t47.*t169;
t180 = t52.*t170;
t181 = t117+t143;
t182 = t122+t138;
t183 = t119+t146;
t184 = t124+t140;
t185 = t121+t149;
t186 = t126+t142;
t188 = t137+t145;
t190 = t139+t148;
t192 = t141+t151;
t207 = -t33.*(t123-t137);
t209 = -t45.*(t125-t139);
t211 = -t50.*(t127-t141);
t222 = -t4.*(t123-t137);
t224 = -t16.*(t125-t139);
t226 = -t21.*(t127-t141);
t250 = t158+t171;
t257 = t162+t172;
t258 = t163+t173;
t259 = t164+t174;
t280 = -t6.*(t153-t172);
t282 = -t18.*(t155-t173);
t283 = -t23.*(t157-t174);
t284 = -t35.*(t153-t172);
t285 = -t47.*(t155-t173);
t286 = -t52.*(t157-t174);
t293 = -t6.*(t95+t4.*(t123-t137));
t295 = -t18.*(t96+t16.*(t125-t139));
t297 = -t23.*(t97+t21.*(t127-t141));
t298 = -t35.*(t95+t4.*(t123-t137));
t299 = -t47.*(t96+t16.*(t125-t139));
t300 = -t52.*(t97+t21.*(t127-t141));
t302 = t6.*(t153-t172);
t303 = t18.*(t155-t173);
t304 = t23.*(t157-t174);
t306 = t6.*(t95+t4.*(t123-t137));
t307 = t18.*(t96+t16.*(t125-t139));
t308 = t23.*(t97+t21.*(t127-t141));
t309 = -t301;
t310 = t172.*1.686712675310875e-2;
t314 = t4.*(t123-t137).*(-1.277082454163948e-1);
t187 = t116+t159;
t189 = t118+t160;
t191 = t120+t161;
t193 = t34.*t181;
t194 = t34.*t182;
t196 = t5.*t181;
t197 = t5.*t182;
t198 = t17.*t183;
t199 = t17.*t184;
t200 = t22.*t185;
t201 = t22.*t186;
t254 = t84+t207;
t255 = t86+t209;
t256 = t88+t211;
t263 = t106+t222;
t264 = t107+t224;
t265 = t108+t226;
t266 = t7.*t250;
t267 = t36.*t250;
t318 = t175+t284;
t319 = t176+t285;
t320 = t177+t286;
t321 = t178+t302;
t322 = t179+t303;
t323 = t180+t304;
t206 = t33.*t187;
t208 = t45.*t189;
t210 = t50.*t191;
t221 = t4.*t187;
t223 = t16.*t189;
t225 = t21.*t191;
t271 = t5.*t254;
t275 = t34.*t254;
t276 = t46.*t255;
t277 = t51.*t256;
t281 = -t266;
t311 = t196.*1.686712675310875e-2;
t312 = t197.*1.686712675310875e-2;
t317 = t267.*2.655364129333133e-1;
t324 = t7.*t318;
t325 = t36.*t318;
t328 = t8.*t321;
t329 = t37.*t321;
t251 = t85+t221;
t252 = t87+t223;
t253 = t89+t225;
t260 = t99+t206;
t261 = t100+t208;
t262 = t101+t210;
t287 = -t5.*(t74-t206);
t291 = -t34.*(t74-t206);
t292 = -t271;
t294 = -t46.*(t75-t208);
t296 = -t51.*(t76-t210);
t305 = t5.*(t74-t206);
t313 = t221.*1.277082454163948e-1;
t315 = -t311;
t316 = -t312;
t326 = t275.*1.686712675310875e-2;
t330 = t34.*(t74-t206).*(-1.686712675310875e-2);
t331 = t34.*(t74-t206).*1.686712675310875e-2;
t332 = t196+t275;
t333 = t198+t276;
t334 = t200+t277;
t356 = t324.*2.655364129333133e-1;
t357 = t329.*2.01340796619765e-2;
t358 = t328.*3.793377327618762e-2;
t361 = t267+t324;
t362 = t281+t325;
t367 = -t8.*(t266-t325);
t370 = -t37.*(t266-t325);
t371 = t8.*(t266-t325);
t396 = t37.*(t266-t325).*(-3.793377327618762e-2);
t268 = t6.*t251;
t269 = t18.*t252;
t270 = t23.*t253;
t272 = t35.*t251;
t273 = t47.*t252;
t274 = t52.*t253;
t327 = -t326;
t335 = t197+t291;
t336 = t199+t294;
t337 = t201+t296;
t338 = t193+t292;
t339 = t6.*t332;
t340 = t18.*t333;
t341 = t23.*t334;
t342 = t35.*t332;
t343 = t47.*t333;
t344 = t52.*t334;
t345 = t194+t305;
t363 = t9.*t361;
t364 = t12.*t361;
t365 = t38.*t361;
t366 = t41.*t361;
t395 = t371.*(-2.01340796619765e-2);
t405 = t328+t370;
t406 = t329+t371;
t288 = -t268;
t289 = -t269;
t290 = -t270;
t346 = t6.*t335;
t347 = t18.*t336;
t348 = t23.*t337;
t349 = t35.*t335;
t350 = t7.*t338;
t351 = t47.*t336;
t352 = t52.*t337;
t353 = t36.*t338;
t354 = t7.*t345;
t355 = t36.*t345;
t368 = -t363;
t369 = -t364;
t372 = t298+t339;
t373 = t299+t340;
t374 = t300+t341;
t378 = t306+t342;
t379 = t307+t343;
t380 = t308+t344;
t407 = t9.*t405;
t408 = t12.*t405;
t409 = t38.*t405;
t410 = t41.*t405;
t411 = t10.*t406;
t412 = t13.*t406;
t413 = t39.*t406;
t414 = t42.*t406;
t359 = t353.*2.655364129333133e-1;
t360 = t355.*2.655364129333133e-1;
t375 = t272+t346;
t376 = t273+t347;
t377 = t274+t348;
t381 = t288+t349;
t382 = t289+t351;
t383 = t290+t352;
t386 = t8.*t372;
t387 = t37.*t372;
t388 = -t7.*(t268-t349);
t389 = -t36.*(t268-t349);
t390 = t7.*t378;
t392 = t36.*t378;
t393 = t7.*(t268-t349);
t415 = -t413;
t416 = -t414;
t442 = t365+t407;
t443 = t366+t408;
t444 = t368+t409;
t445 = t369+t410;
t450 = -t10.*(t363-t409);
t451 = -t13.*(t364-t410);
t452 = -t39.*(t363-t409);
t453 = -t42.*(t364-t410);
t384 = t8.*t375;
t385 = t37.*t375;
t391 = -t386;
t394 = -t390;
t398 = t387.*2.01340796619765e-2;
t399 = t393.*(-2.655364129333133e-1);
t400 = t390.*2.655364129333133e-1;
t402 = t386.*3.793377327618762e-2;
t403 = t393.*2.655364129333133e-1;
t417 = t350+t392;
t418 = t354+t389;
t422 = t355+t393;
t446 = t11.*t442;
t447 = t14.*t443;
t448 = t40.*t442;
t449 = t43.*t443;
t480 = t411+t452;
t481 = t412+t453;
t482 = t415+t450;
t483 = t416+t451;
t397 = t385.*2.01340796619765e-2;
t401 = t384.*3.793377327618762e-2;
t404 = -t400;
t419 = t353+t394;
t420 = t8.*t417;
t421 = t37.*t417;
t427 = t8.*t418;
t428 = t37.*t418;
t429 = t9.*t422;
t430 = t12.*t422;
t432 = t38.*t422;
t433 = t41.*t422;
t484 = t11.*t480;
t485 = t14.*t481;
t486 = t40.*t480;
t487 = t43.*t481;
t423 = t9.*t419;
t424 = t12.*t419;
t425 = t38.*t419;
t426 = t41.*t419;
t431 = -t428;
t434 = -t429;
t435 = -t430;
t436 = t420.*2.01340796619765e-2;
t437 = t421.*3.793377327618762e-2;
t438 = t427.*2.01340796619765e-2;
t440 = t428.*3.793377327618762e-2;
t454 = t387+t420;
t455 = t385+t427;
t456 = t391+t421;
t466 = -t9.*(t386-t421);
t468 = -t12.*(t386-t421);
t470 = -t38.*(t386-t421);
t471 = -t41.*(t386-t421);
t478 = t9.*(t386-t421);
t479 = t12.*(t386-t421);
t488 = -t486;
t489 = -t487;
t516 = t448+t484;
t517 = t449+t485;
t439 = -t437;
t441 = -t440;
t457 = t10.*t454;
t458 = t13.*t454;
t459 = t39.*t454;
t460 = t42.*t454;
t461 = t384+t431;
t462 = t10.*t455;
t463 = t13.*t455;
t464 = t39.*t455;
t465 = t42.*t455;
t490 = t423+t470;
t491 = t424+t471;
t492 = t425+t478;
t493 = t426+t479;
t518 = t446+t488;
t519 = t447+t489;
t520 = t25.*t516;
t521 = t27.*t517;
t467 = -t457;
t469 = -t458;
t472 = t9.*t461;
t473 = t12.*t461;
t474 = t38.*t461;
t475 = -t464;
t476 = t41.*t461;
t477 = -t465;
t494 = t10.*t490;
t495 = t13.*t491;
t496 = t39.*t490;
t497 = t42.*t491;
t500 = t11.*t492;
t501 = t14.*t493;
t502 = t40.*t492;
t503 = t43.*t493;
t522 = t54.*t518;
t523 = t56.*t519;
t498 = t432+t472;
t499 = t433+t473;
t504 = t434+t474;
t505 = t435+t476;
t510 = -t502;
t511 = -t503;
t512 = -t10.*(t429-t474);
t513 = -t13.*(t430-t476);
t514 = -t39.*(t429-t474);
t515 = -t42.*(t430-t476);
t524 = t459+t494;
t525 = t460+t495;
t526 = t467+t496;
t527 = t469+t497;
t528 = -t11.*(t457-t496);
t529 = -t14.*(t458-t497);
t530 = -t40.*(t457-t496);
t531 = -t43.*(t458-t497);
t548 = -t25.*(t502+t11.*(t457-t496));
t549 = -t27.*(t503+t14.*(t458-t497));
t560 = t520+t522;
t561 = t521+t523;
t506 = t11.*t498;
t507 = t14.*t499;
t508 = t40.*t498;
t509 = t43.*t499;
t532 = t462+t514;
t533 = t463+t515;
t534 = t475+t512;
t535 = t477+t513;
t542 = t500+t530;
t543 = t501+t531;
t544 = t510+t528;
t545 = t511+t529;
t536 = t11.*t532;
t537 = t14.*t533;
t538 = t40.*t532;
t539 = t43.*t533;
t546 = t54.*t542;
t547 = t56.*t543;
t540 = -t538;
t541 = -t539;
t550 = -t546;
t551 = -t547;
t552 = t508+t536;
t553 = t509+t537;
t554 = t506+t540;
t555 = t507+t541;
t556 = t25.*t552;
t557 = t27.*t553;
t562 = t548+t550;
t563 = t549+t551;
t558 = t54.*t554;
t559 = t56.*t555;
t564 = t556+t558;
t565 = t557+t559;
Rcut = reshape([t30,t58.*t59,-t29.*t59,0.0,t29,t58,t59,-t30.*t58,t29.*t30,t24.*t323+t53.*t320,t24.*t374-t53.*t380,t24.*t377+t53.*(t270-t352),t24.*t320-t53.*t323,-t24.*t380-t53.*t374,-t53.*t377+t24.*(t270-t352),t51.*t130+t22.*t167,t51.*t185-t22.*t256,t51.*t186+t22.*(t76-t210),t19.*t322+t48.*t319,t19.*t373-t48.*t379,t19.*t376+t48.*(t269-t351),t19.*t319-t48.*t322,-t19.*t379-t48.*t373,-t48.*t376+t19.*(t269-t351),t46.*t129+t17.*t166,t46.*t183-t17.*t255,t46.*t184+t17.*(t75-t208),t321,t372,t375,t361,t419,t422,t266-t325,t417,t418,t27.*t519-t56.*t517,-t56.*(t503+t14.*(t458-t497))+t27.*t543,t27.*t555-t56.*t553,t28.*(t414+t13.*(t364-t410))+t57.*t561,t28.*t525+t57.*(t547+t27.*(t503+t14.*(t458-t497))),t28.*(t465+t13.*(t430-t476))+t57.*t565,-t57.*(t414+t13.*(t364-t410))+t28.*t561,-t57.*t525+t28.*(t547+t27.*(t503+t14.*(t458-t497))),-t57.*(t465+t13.*(t430-t476))+t28.*t565,t25.*t518-t54.*t516,-t54.*(t502+t11.*(t457-t496))+t25.*t542,t25.*t554-t54.*t552,t26.*(t413+t10.*(t363-t409))+t55.*t560,t26.*t524+t55.*(t546+t25.*(t502+t11.*(t457-t496))),t26.*(t464+t10.*(t429-t474))+t55.*t564,-t55.*(t413+t10.*(t363-t409))+t26.*t560,-t55.*t524+t26.*(t546+t25.*(t502+t11.*(t457-t496))),-t55.*(t464+t10.*(t429-t474))+t26.*t564],[3,3,6]);
if nargout > 1
    pcut = reshape([q46,q47,q48,q46-t69.*4.538309500696427e-1+t133.*4.538309500696427e-1-t157.*2.493576648734301e-2+t174.*2.493576648734301e-2-t177.*4.245306301873741e-1-t195+t214+t234-t51.*t130.*9.974306594937201e-3-t22.*t167.*9.974306594937201e-3+t52.*(t157-t174).*4.245306301873741e-1,q47+t97.*4.538309500696427e-1-t200.*2.493576648734301e-2+t219+t231+t235+t241+t246-t277.*2.493576648734301e-2+t308.*4.245306301873741e-1+t344.*4.245306301873741e-1-t51.*t185.*9.974306594937201e-3+t22.*t256.*9.974306594937201e-3+t21.*(t127-t141).*4.538309500696427e-1,q48-t89.*4.538309500696427e-1-t201.*2.493576648734301e-2-t217-t225.*4.538309500696427e-1+t233+t237+t243+t245-t270.*4.245306301873741e-1+t352.*4.245306301873741e-1-t51.*t186.*9.974306594937201e-3-t22.*(t76-t210).*9.974306594937201e-3+t51.*(t76-t210).*2.493576648734301e-2,q46-t68.*4.510689935490981e-1+t132.*4.510689935490981e-1-t155.*2.478401063456583e-2+t173.*2.478401063456583e-2-t176.*4.246654095945845e-1+t195+t214+t234+t46.*t129.*9.913604253826329e-3+t17.*t166.*9.913604253826329e-3+t47.*(t155-t173).*4.246654095945845e-1,q47+t96.*4.510689935490981e-1-t198.*2.478401063456583e-2-t219+t231+t235+t241+t246-t276.*2.478401063456583e-2+t307.*4.246654095945845e-1+t343.*4.246654095945845e-1+t46.*t183.*9.913604253826329e-3-t17.*t255.*9.913604253826329e-3+t16.*(t125-t139).*4.510689935490981e-1,q48-t87.*4.510689935490981e-1-t199.*2.478401063456583e-2+t217-t223.*4.510689935490981e-1+t233+t237+t243+t245-t269.*4.246654095945845e-1+t351.*4.246654095945845e-1+t46.*t184.*9.913604253826329e-3+t17.*(t75-t208).*9.913604253826329e-3+t46.*(t75-t208).*2.478401063456583e-2,q46+t205+t215+t236+t279+t309+t310,q47+t202+t238+t239+t242+t248+t314+t315+t327,q48+t203+t227+t228+t230+t249+t313+t316+t331,q46+t205+t215+t236+t279+t309+t310+t317+t356-t357+t358+t395+t396-t414.*1.371451803080725e-2-t447.*4.065722366579737e-2-t449.*1.674921989294332e-1-t485.*1.674921989294332e-1+t487.*4.065722366579737e-2-t13.*(t364-t410).*1.371451803080725e-2,q47+t202+t238+t239+t242+t248+t314+t315+t327+t359-t398+t402+t404-t436+t439-t460.*1.371451803080725e-2-t495.*1.371451803080725e-2-t501.*4.065722366579737e-2-t503.*1.674921989294332e-1-t14.*(t458-t497).*1.674921989294332e-1+t43.*(t458-t497).*4.065722366579737e-2,q48+t203+t227+t228+t230+t249+t313+t316+t331+t360-t397+t401+t403-t438+t441-t465.*1.371451803080725e-2-t507.*4.065722366579737e-2-t509.*1.674921989294332e-1-t537.*1.674921989294332e-1+t539.*4.065722366579737e-2-t13.*(t430-t476).*1.371451803080725e-2,q46+t205+t215+t236+t279+t309+t310+t317+t356+t357+t358+t371.*2.01340796619765e-2+t396-t413.*1.37145180308193e-2-t446.*4.065722366583309e-2+t448.*1.674921989295803e-1+t484.*1.674921989295803e-1+t486.*4.065722366583309e-2-t10.*(t363-t409).*1.37145180308193e-2,q47+t202+t238+t239+t242+t248+t314+t315+t327+t359+t398+t402+t404+t436+t439-t459.*1.37145180308193e-2-t494.*1.37145180308193e-2-t500.*4.065722366583309e-2+t502.*1.674921989295803e-1+t11.*(t457-t496).*1.674921989295803e-1+t40.*(t457-t496).*4.065722366583309e-2,q48+t203+t227+t228+t230+t249+t313+t316+t331+t360+t397+t401+t403+t438+t441-t464.*1.37145180308193e-2-t506.*4.065722366583309e-2+t508.*1.674921989295803e-1+t536.*1.674921989295803e-1+t538.*4.065722366583309e-2-t10.*(t429-t474).*1.37145180308193e-2],[3,1,6]);
end
