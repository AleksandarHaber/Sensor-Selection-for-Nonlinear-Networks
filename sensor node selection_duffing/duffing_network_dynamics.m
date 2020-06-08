function dx = duffing_network_dynamics(in1)
%DUFFING_NETWORK_DYNAMICS
%    DX = DUFFING_NETWORK_DYNAMICS(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    26-May-2020 00:28:37

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
x5 = in1(5,:);
x6 = in1(6,:);
x7 = in1(7,:);
x8 = in1(8,:);
x9 = in1(9,:);
x10 = in1(10,:);
x11 = in1(11,:);
x12 = in1(12,:);
x13 = in1(13,:);
x14 = in1(14,:);
x15 = in1(15,:);
x16 = in1(16,:);
x17 = in1(17,:);
x18 = in1(18,:);
x19 = in1(19,:);
x20 = in1(20,:);
t2 = -x7;
t3 = -x9;
t4 = -x11;
t5 = -x13;
t6 = -x15;
t7 = -x17;
t8 = -x19;
t9 = t2+x1;
t10 = t3+x1;
t11 = t4+x1;
t12 = t4+x3;
t13 = t6+x1;
t14 = t5+x3;
t15 = t4+x5;
t16 = t8+x1;
t17 = t7+x3;
t18 = t6+x5;
t19 = t7+x5;
t20 = t6+x7;
t21 = t8+x7;
t22 = t6+x11;
t23 = t7+x11;
t24 = t7+x15;
t25 = t8+x15;
t26 = t8+x17;
t27 = t9.^3;
t28 = t10.^3;
t29 = t11.^3;
t30 = t12.^3;
t31 = t13.^3;
t32 = t14.^3;
t33 = t15.^3;
t34 = t16.^3;
t35 = t17.^3;
t36 = t18.^3;
t37 = t19.^3;
t38 = t20.^3;
t39 = t21.^3;
t40 = t22.^3;
t41 = t23.^3;
t42 = t24.^3;
t43 = t25.^3;
t44 = t26.^3;
t45 = t30.*1.235118179871098;
t46 = t39.*1.513212901795988;
t47 = t40.*1.144062849001832;
t48 = t38.*1.610993672302896;
t49 = t36.*1.614230825055887;
t50 = t43.*1.451686750164856;
t51 = t33.*1.32721725212743;
t52 = t28.*1.37038632768104;
t53 = t37.*1.790621028960558;
t54 = t32.*1.786247558529848;
t55 = t27.*1.666649058278422;
t56 = t42.*1.501345910368623;
t57 = t41.*1.478997580336003;
t58 = t31.*1.775515307119208;
t59 = t29.*1.169857792669682;
t60 = t35.*1.45511409891732;
t61 = t34.*1.511434041117034;
t62 = t44.*1.413142738362204;
dx = [x2;t52+t55+t58+t59+t61-x1.*8.992598654384184e+1-x2.*8.064393068563163+x7.*1.387465506109875e+1+x8.*1.489553171500063+x9.*1.526542780364521e+1+x10.*1.079165573117607+x11.*1.704916511853971e+1+x12.*1.460738218389512+x15.*1.285019135146112e+1+x16.*1.485856015669526+x19.*1.337387661604053e+1+x20.*1.511341223646903+x1.^3.*1.780252068321138;x4;t45+t54+t60-x3.*5.669481320421385e+1-x4.*6.815384189614928+x11.*1.216902571990871e+1+x12.*1.711157751533427+x13.*1.579246665643131e+1+x14.*1.640930153327262+x17.*1.619049903815853e+1+x18.*1.809538936085679+x3.^3.*1.353158571222071;x6;t49+t51+t53-x5.*5.676159260145808e+1-x6.*6.157990963929766+x11.*1.257407993768596e+1+x12.*1.791420601593211+x15.*1.467788061290674e+1+x16.*1.437464715641116+x17.*1.365699113933814e+1+x18.*1.429083211105239+x5.^3.*1.744692807074156;x8;t46+t48-t55+x1.*1.387465506109875e+1+x2.*1.489553171500063-x7.*6.174184432062248e+1-x8.*5.871713841785601+x15.*1.366596778790098e+1+x16.*1.848379098029034+x19.*1.640954917060264e+1+x20.*1.505107419792398+x7.^3.*1.435858588580919;x10;-t52+x1.*1.526542780364521e+1+x2.*1.079165573117607-x9.*2.692191509864302e+1-x10.*2.600815415581891+x9.^3.*1.350727103576883;x12;-t45+t47-t51+t57-t59+x1.*1.704916511853971e+1+x2.*1.460738218389512+x3.*1.216902571990871e+1+x4.*1.711157751533427+x5.*1.257407993768596e+1+x6.*1.791420601593211-x11.*8.648216802992773e+1-x12.*9.136032065342285+x15.*1.566101637667775e+1+x16.*1.465765070165355+x17.*1.077071110222013e+1+x18.*1.253152714933861+x11.^3.*1.194764289567049;x14;-t54+x3.*1.579246665643131e+1+x4.*1.640930153327262-x13.*2.663682511154041e+1-x14.*3.167805983835558+x13.^3.*1.438869973126103;x16;-t47-t48-t49+t50+t56-t58+x1.*1.285019135146112e+1+x2.*1.485856015669526+x5.*1.467788061290674e+1+x6.*1.437464715641116+x7.*1.366596778790098e+1+x8.*1.848379098029034+x11.*1.566101637667775e+1+x12.*1.465765070165355-x15.*1.012430409091667e+2-x16.*1.038554951974995e+1+x17.*1.336589568335749e+1+x18.*1.716212311146479+x19.*1.522504322320695e+1+x20.*1.233753906555465+x15.^3.*1.318778301925882;x18;-t53-t56-t57-t60+t62+x3.*1.619049903815853e+1+x4.*1.809538936085679+x5.*1.365699113933814e+1+x6.*1.429083211105239+x11.*1.077071110222013e+1+x12.*1.253152714933861+x15.*1.336589568335749e+1+x16.*1.716212311146479-x17.*7.886917137982404e+1-x18.*9.616797708203529+x19.*1.304599653392559e+1+x20.*1.991066430615609+x17.^3.*1.237283579771522;x20;-t46-t50-t61-t62+x1.*1.337387661604053e+1+x2.*1.511341223646903+x7.*1.640954917060264e+1+x8.*1.505107419792398+x15.*1.522504322320695e+1+x16.*1.233753906555465+x17.*1.304599653392559e+1+x18.*1.991066430615609-x19.*6.916649309671358e+1-x20.*7.412390046966806+x19.^3.*1.987982003161633];