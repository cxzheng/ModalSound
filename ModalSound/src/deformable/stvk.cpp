#include "stvk.h"
#include "tensor.hpp"

/*!
 * the Cauchy stress tensor is:
 * FSF^T / det[F]
 */
void StVKMaterial::cauchy_stress_tensor(
        const Matrix3<REAL>& F, Matrix3<REAL>& out) const
{
    second_piola_kirchhoff(F, out); // out = S

    out = F * out * F.transpose();
    out /= F.det();
}

void StVKMaterial::first_piola_kirchhoff(
        const Matrix3<REAL>& F, Matrix3<REAL>& out) const
{
    second_piola_kirchhoff(F, out);
    out = F * out;
}

void StVKMaterial::second_piola_kirchhoff(
        const Matrix3<REAL>& F, Matrix3<REAL>& out) const
{
    Matrix3<REAL> C, E;
    right_cauchy_green_deformation_tensor(F, C);
    green_stain_tensor(C, E);

    out = Matrix3<REAL>::identity(lambda_ * E.trace()) + 
        (REAL)2 * mu_ * E;
}

REAL StVKMaterial::strain_energy_density(const Matrix3<REAL>& F) const
{
    Matrix3<REAL> C, E;
    right_cauchy_green_deformation_tensor(F, C);
    green_stain_tensor(C, E);

    const REAL tr = E.trace();
    return 0.5*lambda_*M_SQR(tr) + mu_*E.colon(E); 
}

/*
 * Stiffness matrix: generated code from MATLAB
 */
void StVKMaterial::stiffness_matrix(const Tet<REAL>& tet,
        Matrix<REAL>& out) const
{
    const Vector3<REAL>* b = tet.weighted_normals();
    const Tet<REAL>::TVtx& n0 = tet[0];
    const Tet<REAL>::TVtx& n1 = tet[1];
    const Tet<REAL>::TVtx& n2 = tet[2];
    const Tet<REAL>::TVtx& n3 = tet[3];
    const Matrix3<REAL>& iDm = tet.inverse_Dm();

    const REAL t1 = -iDm(0,0)-iDm(1,0)-iDm(2,0);
    const REAL t2 = n1.x-n0.x;
    const REAL t4 = n2.x-n0.x;
    const REAL t6 = n3.x-n0.x;
    const REAL t8 = t2*iDm(0,0)+t4*iDm(1,0)+t6*iDm(2,0);
    const REAL t9 = t8*t8;
    const REAL t11 = n1.y-n0.y;
    const REAL t13 = n2.y-n0.y;
    const REAL t15 = n3.y-n0.y;
    const REAL t17 = t11*iDm(0,0)+t13*iDm(1,0)+t15*iDm(2,0);
    const REAL t18 = t17*t17;
    const REAL t20 = n1.z-n0.z;
    const REAL t22 = n2.z-n0.z;
    const REAL t24 = n3.z-n0.z;
    const REAL t26 = t20*iDm(0,0)+t22*iDm(1,0)+t24*iDm(2,0);
    const REAL t27 = t26*t26;
    const REAL t32 = t2*iDm(0,1)+t4*iDm(1,1)+t6*iDm(2,1);
    const REAL t33 = t32*t32;
    const REAL t38 = t11*iDm(0,1)+t13*iDm(1,1)+t15*iDm(2,1);
    const REAL t39 = t38*t38;
    const REAL t44 = t20*iDm(0,1)+t22*iDm(1,1)+t24*iDm(2,1);
    const REAL t45 = t44*t44;
    const REAL t50 = t2*iDm(0,2)+t4*iDm(1,2)+t6*iDm(2,2);
    const REAL t51 = t50*t50;
    const REAL t56 = t11*iDm(0,2)+t13*iDm(1,2)+t15*iDm(2,2);
    const REAL t57 = t56*t56;
    const REAL t62 = t20*iDm(0,2)+t22*iDm(1,2)+t24*iDm(2,2);
    const REAL t63 = t62*t62;
    const REAL t66 = lambda_*(t9/2.0+t18/2.0+t27/2.0-3.0/2.0+t33/2.0+t39/2.0+t45/2.0+t51/2.0+t57/2.0+t63/2.0);
    const REAL t70 = t66+mu_*(t9+t18+t27-1.0);
    const REAL t71 = t1*t70;
    const REAL t73 = -iDm(0,1)-iDm(1,1)-iDm(2,1);
    const REAL t75 = -iDm(0,2)-iDm(1,2)-iDm(2,2);
    const REAL t78 = lambda_*(t8*t1+t32*t73+t50*t75);
    const REAL t79 = mu_*t8;
    const REAL t84 = t73*mu_;
    const REAL t88 = t8*t32+t17*t38+t26*t44;
    const REAL t90 = t84*t88;
    const REAL t91 = t32*mu_;
    const REAL t94 = t1*t32+t8*t73;
    const REAL t97 = t75*mu_;
    const REAL t101 = t8*t50+t17*t56+t26*t62;
    const REAL t103 = t97*t101;
    const REAL t104 = t50*mu_;
    const REAL t107 = t1*t50+t8*t75;
    const REAL t112 = t1*mu_;
    const REAL t114 = t112*t88;
    const REAL t120 = t66+mu_*(t33+t39+t45-1.0);
    const REAL t121 = t73*t120;
    const REAL t129 = t32*t50+t38*t56+t44*t62;
    const REAL t131 = t97*t129;
    const REAL t134 = t73*t50+t32*t75;
    const REAL t140 = t112*t101;
    const REAL t144 = t84*t129;
    const REAL t150 = t66+mu_*(t51+t57+t63-1.0);
    const REAL t151 = t75*t150;
    const REAL t163 = lambda_*(t17*t1+t38*t73+t56*t75);
    const REAL t164 = mu_*t17;
    const REAL t167 = t163+2.0*t164*t1;
    const REAL t171 = t1*t38+t17*t73;
    const REAL t176 = t1*t56+t17*t75;
    const REAL t183 = mu_*t38;
    const REAL t186 = t163+2.0*t183*t73;
    const REAL t190 = t73*t56+t38*t75;
    const REAL t199 = mu_*t56;
    const REAL t202 = t163+2.0*t199*t75;
    const REAL t211 = lambda_*(t26*t1+t44*t73+t62*t75);
    const REAL t212 = mu_*t26;
    const REAL t215 = t211+2.0*t212*t1;
    const REAL t219 = t1*t44+t26*t73;
    const REAL t224 = t1*t62+t26*t75;
    const REAL t231 = mu_*t44;
    const REAL t234 = t211+2.0*t231*t73;
    const REAL t238 = t73*t62+t44*t75;
    const REAL t247 = mu_*t62;
    const REAL t250 = t211+2.0*t247*t75;
    const REAL t255 = iDm(0,0)*t70;
    const REAL t260 = lambda_*(t8*iDm(0,0)+t32*iDm(0,1)+t50*iDm(0,2));
    const REAL t263 = t260+2.0*t79*iDm(0,0);
    const REAL t265 = iDm(0,1)*mu_;
    const REAL t267 = t265*t88;
    const REAL t270 = iDm(0,0)*t32+t8*iDm(0,1);
    const REAL t273 = iDm(0,2)*mu_;
    const REAL t275 = t273*t101;
    const REAL t278 = iDm(0,0)*t50+t8*iDm(0,2);
    const REAL t281 = t255+t8*t263+t267+t91*t270+t275+t104*t278;
    const REAL t283 = iDm(0,0)*mu_;
    const REAL t285 = t283*t88;
    const REAL t288 = iDm(0,1)*t120;
    const REAL t291 = t260+2.0*t91*iDm(0,1);
    const REAL t294 = t273*t129;
    const REAL t297 = iDm(0,1)*t50+t32*iDm(0,2);
    const REAL t300 = t285+t79*t270+t288+t32*t291+t294+t104*t297;
    const REAL t303 = t283*t101;
    const REAL t307 = t265*t129;
    const REAL t310 = iDm(0,2)*t150;
    const REAL t313 = t260+2.0*t104*iDm(0,2);
    const REAL t315 = t303+t79*t278+t307+t91*t297+t310+t50*t313;
    const REAL t322 = lambda_*(t17*iDm(0,0)+t38*iDm(0,1)+t56*iDm(0,2));
    const REAL t325 = t322+2.0*t164*iDm(0,0);
    const REAL t329 = iDm(0,0)*t38+t17*iDm(0,1);
    const REAL t334 = iDm(0,0)*t56+t17*iDm(0,2);
    const REAL t337 = t8*t325+t91*t329+t104*t334;
    const REAL t343 = t322+2.0*t183*iDm(0,1);
    const REAL t347 = iDm(0,1)*t56+t38*iDm(0,2);
    const REAL t350 = t79*t329+t32*t343+t104*t347;
    const REAL t358 = t322+2.0*t199*iDm(0,2);
    const REAL t360 = t79*t334+t91*t347+t50*t358;
    const REAL t367 = lambda_*(t26*iDm(0,0)+t44*iDm(0,1)+t62*iDm(0,2));
    const REAL t370 = t367+2.0*t212*iDm(0,0);
    const REAL t374 = iDm(0,0)*t44+t26*iDm(0,1);
    const REAL t379 = iDm(0,0)*t62+t26*iDm(0,2);
    const REAL t382 = t8*t370+t91*t374+t104*t379;
    const REAL t388 = t367+2.0*t231*iDm(0,1);
    const REAL t392 = iDm(0,1)*t62+t44*iDm(0,2);
    const REAL t395 = t79*t374+t32*t388+t104*t392;
    const REAL t403 = t367+2.0*t247*iDm(0,2);
    const REAL t405 = t79*t379+t91*t392+t50*t403;
    const REAL t408 = iDm(1,0)*t70;
    const REAL t413 = lambda_*(t8*iDm(1,0)+t32*iDm(1,1)+t50*iDm(1,2));
    const REAL t416 = t413+2.0*t79*iDm(1,0);
    const REAL t418 = iDm(1,1)*mu_;
    const REAL t420 = t418*t88;
    const REAL t423 = iDm(1,0)*t32+t8*iDm(1,1);
    const REAL t426 = iDm(1,2)*mu_;
    const REAL t428 = t426*t101;
    const REAL t431 = iDm(1,0)*t50+t8*iDm(1,2);
    const REAL t434 = t408+t8*t416+t420+t91*t423+t428+t104*t431;
    const REAL t436 = iDm(1,0)*mu_;
    const REAL t438 = t436*t88;
    const REAL t441 = iDm(1,1)*t120;
    const REAL t444 = t413+2.0*t91*iDm(1,1);
    const REAL t447 = t426*t129;
    const REAL t450 = iDm(1,1)*t50+t32*iDm(1,2);
    const REAL t453 = t438+t79*t423+t441+t32*t444+t447+t104*t450;
    const REAL t456 = t436*t101;
    const REAL t460 = t418*t129;
    const REAL t463 = iDm(1,2)*t150;
    const REAL t466 = t413+2.0*t104*iDm(1,2);
    const REAL t468 = t456+t79*t431+t460+t91*t450+t463+t50*t466;
    const REAL t475 = lambda_*(t17*iDm(1,0)+t38*iDm(1,1)+t56*iDm(1,2));
    const REAL t478 = t475+2.0*t164*iDm(1,0);
    const REAL t482 = iDm(1,0)*t38+t17*iDm(1,1);
    const REAL t487 = iDm(1,0)*t56+t17*iDm(1,2);
    const REAL t490 = t8*t478+t91*t482+t104*t487;
    const REAL t496 = t475+2.0*t183*iDm(1,1);
    const REAL t500 = iDm(1,1)*t56+t38*iDm(1,2);
    const REAL t503 = t79*t482+t32*t496+t104*t500;
    const REAL t511 = t475+2.0*t199*iDm(1,2);
    const REAL t513 = t79*t487+t91*t500+t50*t511;
    const REAL t520 = lambda_*(t26*iDm(1,0)+t44*iDm(1,1)+t62*iDm(1,2));
    const REAL t523 = t520+2.0*t212*iDm(1,0);
    const REAL t527 = iDm(1,0)*t44+t26*iDm(1,1);
    const REAL t532 = iDm(1,0)*t62+t26*iDm(1,2);
    const REAL t535 = t8*t523+t91*t527+t104*t532;
    const REAL t541 = t520+2.0*t231*iDm(1,1);
    const REAL t545 = iDm(1,1)*t62+t44*iDm(1,2);
    const REAL t548 = t79*t527+t32*t541+t104*t545;
    const REAL t556 = t520+2.0*t247*iDm(1,2);
    const REAL t558 = t79*t532+t91*t545+t50*t556;
    const REAL t561 = iDm(2,0)*t70;
    const REAL t566 = lambda_*(t8*iDm(2,0)+t32*iDm(2,1)+t50*iDm(2,2));
    const REAL t569 = t566+2.0*t79*iDm(2,0);
    const REAL t571 = iDm(2,1)*mu_;
    const REAL t573 = t571*t88;
    const REAL t576 = iDm(2,0)*t32+t8*iDm(2,1);
    const REAL t579 = iDm(2,2)*mu_;
    const REAL t581 = t579*t101;
    const REAL t584 = iDm(2,0)*t50+t8*iDm(2,2);
    const REAL t587 = t561+t8*t569+t573+t91*t576+t581+t104*t584;
    const REAL t589 = iDm(2,0)*mu_;
    const REAL t591 = t589*t88;
    const REAL t594 = iDm(2,1)*t120;
    const REAL t597 = t566+2.0*t91*iDm(2,1);
    const REAL t600 = t579*t129;
    const REAL t603 = iDm(2,1)*t50+t32*iDm(2,2);
    const REAL t606 = t591+t79*t576+t594+t32*t597+t600+t104*t603;
    const REAL t609 = t589*t101;
    const REAL t613 = t571*t129;
    const REAL t616 = iDm(2,2)*t150;
    const REAL t619 = t566+2.0*t104*iDm(2,2);
    const REAL t621 = t609+t79*t584+t613+t91*t603+t616+t50*t619;
    const REAL t628 = lambda_*(t17*iDm(2,0)+t38*iDm(2,1)+t56*iDm(2,2));
    const REAL t631 = t628+2.0*t164*iDm(2,0);
    const REAL t635 = iDm(2,0)*t38+t17*iDm(2,1);
    const REAL t640 = iDm(2,0)*t56+t17*iDm(2,2);
    const REAL t643 = t8*t631+t91*t635+t104*t640;
    const REAL t649 = t628+2.0*t183*iDm(2,1);
    const REAL t653 = iDm(2,1)*t56+t38*iDm(2,2);
    const REAL t656 = t79*t635+t32*t649+t104*t653;
    const REAL t664 = t628+2.0*t199*iDm(2,2);
    const REAL t666 = t79*t640+t91*t653+t50*t664;
    const REAL t673 = lambda_*(t26*iDm(2,0)+t44*iDm(2,1)+t62*iDm(2,2));
    const REAL t676 = t673+2.0*t212*iDm(2,0);
    const REAL t680 = iDm(2,0)*t44+t26*iDm(2,1);
    const REAL t685 = iDm(2,0)*t62+t26*iDm(2,2);
    const REAL t688 = t8*t676+t91*t680+t104*t685;
    const REAL t694 = t673+2.0*t231*iDm(2,1);
    const REAL t698 = iDm(2,1)*t62+t44*iDm(2,2);
    const REAL t701 = t79*t680+t32*t694+t104*t698;
    const REAL t709 = t673+2.0*t247*iDm(2,2);
    const REAL t711 = t79*t685+t91*t698+t50*t709;
    const REAL t785 = t255+t17*t325+t267+t183*t329+t275+t199*t334;
    const REAL t792 = t285+t164*t329+t288+t38*t343+t294+t199*t347;
    const REAL t799 = t303+t164*t334+t307+t183*t347+t310+t56*t358;
    const REAL t807 = t17*t370+t183*t374+t199*t379;
    const REAL t814 = t164*t374+t38*t388+t199*t392;
    const REAL t821 = t164*t379+t183*t392+t56*t403;
    const REAL t829 = t17*t416+t183*t423+t199*t431;
    const REAL t836 = t164*t423+t38*t444+t199*t450;
    const REAL t843 = t164*t431+t183*t450+t56*t466;
    const REAL t851 = t408+t17*t478+t420+t183*t482+t428+t199*t487;
    const REAL t858 = t438+t164*t482+t441+t38*t496+t447+t199*t500;
    const REAL t865 = t456+t164*t487+t460+t183*t500+t463+t56*t511;
    const REAL t873 = t17*t523+t183*t527+t199*t532;
    const REAL t880 = t164*t527+t38*t541+t199*t545;
    const REAL t887 = t164*t532+t183*t545+t56*t556;
    const REAL t895 = t17*t569+t183*t576+t199*t584;
    const REAL t902 = t164*t576+t38*t597+t199*t603;
    const REAL t909 = t164*t584+t183*t603+t56*t619;
    const REAL t917 = t561+t17*t631+t573+t183*t635+t581+t199*t640;
    const REAL t924 = t591+t164*t635+t594+t38*t649+t600+t199*t653;
    const REAL t931 = t609+t164*t640+t613+t183*t653+t616+t56*t664;
    const REAL t939 = t17*t676+t183*t680+t199*t685;
    const REAL t946 = t164*t680+t38*t694+t199*t698;
    const REAL t953 = t164*t685+t183*t698+t56*t709;
    const REAL t1027 = t255+t26*t370+t267+t231*t374+t275+t247*t379;
    const REAL t1034 = t285+t212*t374+t288+t44*t388+t294+t247*t392;
    const REAL t1041 = t303+t212*t379+t307+t231*t392+t310+t62*t403;
    const REAL t1049 = t26*t416+t231*t423+t247*t431;
    const REAL t1056 = t212*t423+t44*t444+t247*t450;
    const REAL t1063 = t212*t431+t231*t450+t62*t466;
    const REAL t1071 = t26*t478+t231*t482+t247*t487;
    const REAL t1078 = t212*t482+t44*t496+t247*t500;
    const REAL t1085 = t212*t487+t231*t500+t62*t511;
    const REAL t1093 = t408+t26*t523+t420+t231*t527+t428+t247*t532;
    const REAL t1100 = t438+t212*t527+t441+t44*t541+t447+t247*t545;
    const REAL t1107 = t456+t212*t532+t460+t231*t545+t463+t62*t556;
    const REAL t1115 = t26*t569+t231*t576+t247*t584;
    const REAL t1122 = t212*t576+t44*t597+t247*t603;
    const REAL t1129 = t212*t584+t231*t603+t62*t619;
    const REAL t1137 = t26*t631+t231*t635+t247*t640;
    const REAL t1144 = t212*t635+t44*t649+t247*t653;
    const REAL t1151 = t212*t640+t231*t653+t62*t664;
    const REAL t1159 = t561+t26*t676+t573+t231*t680+t581+t247*t685;
    const REAL t1166 = t591+t212*t680+t594+t44*t694+t600+t247*t698;
    const REAL t1173 = t609+t212*t685+t613+t231*t698+t616+t62*t709;

    out[0][0] = (t71+t8*(t78+2.0*t79*t1)+t90+t91*t94+t103+t104*t107)*b[0].x+(t114+t79*t94+t121+t32*(t78+2.0*t91*t73)+t131+t104*t134)*b[0].y+(t140+t79*t107+t144+t91*t134+t151+t50*(t78+2.0*t104*t75))*b[0].z;
    out[0][1] = (t8*t167+t91*t171+t104*t176)*b[0].x+(t79*t171+t32*t186+t104*t190)*b[0].y+(t79*t176+t91*t190+t50*t202)*b[0].z;
    out[0][2] = (t8*t215+t91*t219+t104*t224)*b[0].x+(t79*t219+t32*t234+t104*t238)*b[0].y+(t79*t224+t91*t238+t50*t250)*b[0].z;
    out[0][3] = t281*b[0].x+t300*b[0].y+t315*b[0].z;
    out[0][4] = t337*b[0].x+t350*b[0].y+t360*b[0].z;
    out[0][5] = t382*b[0].x+t395*b[0].y+t405*b[0].z;
    out[0][6] = t434*b[0].x+t453*b[0].y+t468*b[0].z;
    out[0][7] = t490*b[0].x+t503*b[0].y+t513*b[0].z;
    out[0][8] = t535*b[0].x+t548*b[0].y+t558*b[0].z;
    out[0][9] = t587*b[0].x+t606*b[0].y+t621*b[0].z;
    out[0][10] = t643*b[0].x+t656*b[0].y+t666*b[0].z;
    out[0][11] = t688*b[0].x+t701*b[0].y+t711*b[0].z;
    out[1][1] = (t71+t17*t167+t90+t183*t171+t103+t199*t176)*b[0].x+(t114+t164*t171+t121+t38*t186+t131+t199*t190)*b[0].y+(t140+t164*t176+t144+t183*t190+t151+t56*t202)*b[0].z;
    out[1][2] = (t17*t215+t183*t219+t199*t224)*b[0].x+(t164*t219+t38*t234+t199*t238)*b[0].y+(t164*t224+t183*t238+t56*t250)*b[0].z;
    out[1][3] = (t17*t263+t183*t270+t199*t278)*b[0].x+(t164*t270+t38*t291+t199*t297)*b[0].y+(t164*t278+t183*t297+t56*t313)*b[0].z;
    out[1][4] = t785*b[0].x+t792*b[0].y+t799*b[0].z;
    out[1][5] = t807*b[0].x+t814*b[0].y+t821*b[0].z;
    out[1][6] = t829*b[0].x+t836*b[0].y+t843*b[0].z;
    out[1][7] = t851*b[0].x+t858*b[0].y+t865*b[0].z;
    out[1][8] = t873*b[0].x+t880*b[0].y+t887*b[0].z;
    out[1][9] = t895*b[0].x+t902*b[0].y+t909*b[0].z;
    out[1][10] = t917*b[0].x+t924*b[0].y+t931*b[0].z;
    out[1][11] = t939*b[0].x+t946*b[0].y+t953*b[0].z;
    out[2][2] = (t71+t26*t215+t90+t231*t219+t103+t247*t224)*b[0].x+(t114+t212*t219+t121+t44*t234+t131+t247*t238)*b[0].y+(t140+t212*t224+t144+t231*t238+t151+t62*t250)*b[0].z;
    out[2][3] = (t26*t263+t231*t270+t247*t278)*b[0].x+(t212*t270+t44*t291+t247*t297)*b[0].y+(t212*t278+t231*t297+t62*t313)*b[0].z;
    out[2][4] = (t26*t325+t231*t329+t247*t334)*b[0].x+(t212*t329+t44*t343+t247*t347)*b[0].y+(t212*t334+t231*t347+t62*t358)*b[0].z;
    out[2][5] = t1027*b[0].x+t1034*b[0].y+t1041*b[0].z;
    out[2][6] = t1049*b[0].x+t1056*b[0].y+t1063*b[0].z;
    out[2][7] = t1071*b[0].x+t1078*b[0].y+t1085*b[0].z;
    out[2][8] = t1093*b[0].x+t1100*b[0].y+t1107*b[0].z;
    out[2][9] = t1115*b[0].x+t1122*b[0].y+t1129*b[0].z;
    out[2][10] = t1137*b[0].x+t1144*b[0].y+t1151*b[0].z;
    out[2][11] = t1159*b[0].x+t1166*b[0].y+t1173*b[0].z;
    out[3][3] = t281*b[1].x+t300*b[1].y+t315*b[1].z;
    out[3][4] = t337*b[1].x+t350*b[1].y+t360*b[1].z;
    out[3][5] = t382*b[1].x+t395*b[1].y+t405*b[1].z;
    out[3][6] = t434*b[1].x+t453*b[1].y+t468*b[1].z;
    out[3][7] = t490*b[1].x+t503*b[1].y+t513*b[1].z;
    out[3][8] = t535*b[1].x+t548*b[1].y+t558*b[1].z;
    out[3][9] = t587*b[1].x+t606*b[1].y+t621*b[1].z;
    out[3][10] = t643*b[1].x+t656*b[1].y+t666*b[1].z;
    out[3][11] = t688*b[1].x+t701*b[1].y+t711*b[1].z;
    out[4][4] = t785*b[1].x+t792*b[1].y+t799*b[1].z;
    out[4][5] = t807*b[1].x+t814*b[1].y+t821*b[1].z;
    out[4][6] = t829*b[1].x+t836*b[1].y+t843*b[1].z;
    out[4][7] = t851*b[1].x+t858*b[1].y+t865*b[1].z;
    out[4][8] = t873*b[1].x+t880*b[1].y+t887*b[1].z;
    out[4][9] = t895*b[1].x+t902*b[1].y+t909*b[1].z;
    out[4][10] = t917*b[1].x+t924*b[1].y+t931*b[1].z;
    out[4][11] = t939*b[1].x+t946*b[1].y+t953*b[1].z;
    out[5][5] = t1027*b[1].x+t1034*b[1].y+t1041*b[1].z;
    out[5][6] = t1049*b[1].x+t1056*b[1].y+t1063*b[1].z;
    out[5][7] = t1071*b[1].x+t1078*b[1].y+t1085*b[1].z;
    out[5][8] = t1093*b[1].x+t1100*b[1].y+t1107*b[1].z;
    out[5][9] = t1115*b[1].x+t1122*b[1].y+t1129*b[1].z;
    out[5][10] = t1137*b[1].x+t1144*b[1].y+t1151*b[1].z;
    out[5][11] = t1159*b[1].x+t1166*b[1].y+t1173*b[1].z;
    out[6][6] = t434*b[2].x+t453*b[2].y+t468*b[2].z;
    out[6][7] = t490*b[2].x+t503*b[2].y+t513*b[2].z;
    out[6][8] = t535*b[2].x+t548*b[2].y+t558*b[2].z;
    out[6][9] = t587*b[2].x+t606*b[2].y+t621*b[2].z;
    out[6][10] = t643*b[2].x+t656*b[2].y+t666*b[2].z;
    out[6][11] = t688*b[2].x+t701*b[2].y+t711*b[2].z;
    out[7][7] = t851*b[2].x+t858*b[2].y+t865*b[2].z;
    out[7][8] = t873*b[2].x+t880*b[2].y+t887*b[2].z;
    out[7][9] = t895*b[2].x+t902*b[2].y+t909*b[2].z;
    out[7][10] = t917*b[2].x+t924*b[2].y+t931*b[2].z;
    out[7][11] = t939*b[2].x+t946*b[2].y+t953*b[2].z;
    out[8][8] = t1093*b[2].x+t1100*b[2].y+t1107*b[2].z;
    out[8][9] = t1115*b[2].x+t1122*b[2].y+t1129*b[2].z;
    out[8][10] = t1137*b[2].x+t1144*b[2].y+t1151*b[2].z;
    out[8][11] = t1159*b[2].x+t1166*b[2].y+t1173*b[2].z;
    out[9][9] = t587*b[3].x+t606*b[3].y+t621*b[3].z;
    out[9][10] = t643*b[3].x+t656*b[3].y+t666*b[3].z;
    out[9][11] = t688*b[3].x+t701*b[3].y+t711*b[3].z;
    out[10][10] = t917*b[3].x+t924*b[3].y+t931*b[3].z;
    out[10][11] = t939*b[3].x+t946*b[3].y+t953*b[3].z;
    out[11][11] = t1159*b[3].x+t1166*b[3].y+t1173*b[3].z;
}
