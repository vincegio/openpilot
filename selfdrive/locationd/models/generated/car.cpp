#include "car.h"

namespace {
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 5.991464547107981;

/******************************************************************************
 *                       Code generated with sympy 1.8                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_312213572642374608) {
   out_312213572642374608[0] = delta_x[0] + nom_x[0];
   out_312213572642374608[1] = delta_x[1] + nom_x[1];
   out_312213572642374608[2] = delta_x[2] + nom_x[2];
   out_312213572642374608[3] = delta_x[3] + nom_x[3];
   out_312213572642374608[4] = delta_x[4] + nom_x[4];
   out_312213572642374608[5] = delta_x[5] + nom_x[5];
   out_312213572642374608[6] = delta_x[6] + nom_x[6];
   out_312213572642374608[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_486699098684977001) {
   out_486699098684977001[0] = -nom_x[0] + true_x[0];
   out_486699098684977001[1] = -nom_x[1] + true_x[1];
   out_486699098684977001[2] = -nom_x[2] + true_x[2];
   out_486699098684977001[3] = -nom_x[3] + true_x[3];
   out_486699098684977001[4] = -nom_x[4] + true_x[4];
   out_486699098684977001[5] = -nom_x[5] + true_x[5];
   out_486699098684977001[6] = -nom_x[6] + true_x[6];
   out_486699098684977001[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_3333653717063534221) {
   out_3333653717063534221[0] = 1.0;
   out_3333653717063534221[1] = 0.0;
   out_3333653717063534221[2] = 0.0;
   out_3333653717063534221[3] = 0.0;
   out_3333653717063534221[4] = 0.0;
   out_3333653717063534221[5] = 0.0;
   out_3333653717063534221[6] = 0.0;
   out_3333653717063534221[7] = 0.0;
   out_3333653717063534221[8] = 0.0;
   out_3333653717063534221[9] = 1.0;
   out_3333653717063534221[10] = 0.0;
   out_3333653717063534221[11] = 0.0;
   out_3333653717063534221[12] = 0.0;
   out_3333653717063534221[13] = 0.0;
   out_3333653717063534221[14] = 0.0;
   out_3333653717063534221[15] = 0.0;
   out_3333653717063534221[16] = 0.0;
   out_3333653717063534221[17] = 0.0;
   out_3333653717063534221[18] = 1.0;
   out_3333653717063534221[19] = 0.0;
   out_3333653717063534221[20] = 0.0;
   out_3333653717063534221[21] = 0.0;
   out_3333653717063534221[22] = 0.0;
   out_3333653717063534221[23] = 0.0;
   out_3333653717063534221[24] = 0.0;
   out_3333653717063534221[25] = 0.0;
   out_3333653717063534221[26] = 0.0;
   out_3333653717063534221[27] = 1.0;
   out_3333653717063534221[28] = 0.0;
   out_3333653717063534221[29] = 0.0;
   out_3333653717063534221[30] = 0.0;
   out_3333653717063534221[31] = 0.0;
   out_3333653717063534221[32] = 0.0;
   out_3333653717063534221[33] = 0.0;
   out_3333653717063534221[34] = 0.0;
   out_3333653717063534221[35] = 0.0;
   out_3333653717063534221[36] = 1.0;
   out_3333653717063534221[37] = 0.0;
   out_3333653717063534221[38] = 0.0;
   out_3333653717063534221[39] = 0.0;
   out_3333653717063534221[40] = 0.0;
   out_3333653717063534221[41] = 0.0;
   out_3333653717063534221[42] = 0.0;
   out_3333653717063534221[43] = 0.0;
   out_3333653717063534221[44] = 0.0;
   out_3333653717063534221[45] = 1.0;
   out_3333653717063534221[46] = 0.0;
   out_3333653717063534221[47] = 0.0;
   out_3333653717063534221[48] = 0.0;
   out_3333653717063534221[49] = 0.0;
   out_3333653717063534221[50] = 0.0;
   out_3333653717063534221[51] = 0.0;
   out_3333653717063534221[52] = 0.0;
   out_3333653717063534221[53] = 0.0;
   out_3333653717063534221[54] = 1.0;
   out_3333653717063534221[55] = 0.0;
   out_3333653717063534221[56] = 0.0;
   out_3333653717063534221[57] = 0.0;
   out_3333653717063534221[58] = 0.0;
   out_3333653717063534221[59] = 0.0;
   out_3333653717063534221[60] = 0.0;
   out_3333653717063534221[61] = 0.0;
   out_3333653717063534221[62] = 0.0;
   out_3333653717063534221[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_1478154795453816204) {
   out_1478154795453816204[0] = state[0];
   out_1478154795453816204[1] = state[1];
   out_1478154795453816204[2] = state[2];
   out_1478154795453816204[3] = state[3];
   out_1478154795453816204[4] = state[4];
   out_1478154795453816204[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_1478154795453816204[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_1478154795453816204[7] = state[7];
}
void F_fun(double *state, double dt, double *out_2267919153160654385) {
   out_2267919153160654385[0] = 1;
   out_2267919153160654385[1] = 0;
   out_2267919153160654385[2] = 0;
   out_2267919153160654385[3] = 0;
   out_2267919153160654385[4] = 0;
   out_2267919153160654385[5] = 0;
   out_2267919153160654385[6] = 0;
   out_2267919153160654385[7] = 0;
   out_2267919153160654385[8] = 0;
   out_2267919153160654385[9] = 1;
   out_2267919153160654385[10] = 0;
   out_2267919153160654385[11] = 0;
   out_2267919153160654385[12] = 0;
   out_2267919153160654385[13] = 0;
   out_2267919153160654385[14] = 0;
   out_2267919153160654385[15] = 0;
   out_2267919153160654385[16] = 0;
   out_2267919153160654385[17] = 0;
   out_2267919153160654385[18] = 1;
   out_2267919153160654385[19] = 0;
   out_2267919153160654385[20] = 0;
   out_2267919153160654385[21] = 0;
   out_2267919153160654385[22] = 0;
   out_2267919153160654385[23] = 0;
   out_2267919153160654385[24] = 0;
   out_2267919153160654385[25] = 0;
   out_2267919153160654385[26] = 0;
   out_2267919153160654385[27] = 1;
   out_2267919153160654385[28] = 0;
   out_2267919153160654385[29] = 0;
   out_2267919153160654385[30] = 0;
   out_2267919153160654385[31] = 0;
   out_2267919153160654385[32] = 0;
   out_2267919153160654385[33] = 0;
   out_2267919153160654385[34] = 0;
   out_2267919153160654385[35] = 0;
   out_2267919153160654385[36] = 1;
   out_2267919153160654385[37] = 0;
   out_2267919153160654385[38] = 0;
   out_2267919153160654385[39] = 0;
   out_2267919153160654385[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_2267919153160654385[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_2267919153160654385[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2267919153160654385[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2267919153160654385[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_2267919153160654385[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_2267919153160654385[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_2267919153160654385[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_2267919153160654385[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_2267919153160654385[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_2267919153160654385[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2267919153160654385[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2267919153160654385[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_2267919153160654385[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_2267919153160654385[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_2267919153160654385[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2267919153160654385[56] = 0;
   out_2267919153160654385[57] = 0;
   out_2267919153160654385[58] = 0;
   out_2267919153160654385[59] = 0;
   out_2267919153160654385[60] = 0;
   out_2267919153160654385[61] = 0;
   out_2267919153160654385[62] = 0;
   out_2267919153160654385[63] = 1;
}
void h_25(double *state, double *unused, double *out_2116182985142631742) {
   out_2116182985142631742[0] = state[6];
}
void H_25(double *state, double *unused, double *out_6245485705572964003) {
   out_6245485705572964003[0] = 0;
   out_6245485705572964003[1] = 0;
   out_6245485705572964003[2] = 0;
   out_6245485705572964003[3] = 0;
   out_6245485705572964003[4] = 0;
   out_6245485705572964003[5] = 0;
   out_6245485705572964003[6] = 1;
   out_6245485705572964003[7] = 0;
}
void h_24(double *state, double *unused, double *out_5094452919565144864) {
   out_5094452919565144864[0] = state[4];
   out_5094452919565144864[1] = state[5];
}
void H_24(double *state, double *unused, double *out_7507874769768704135) {
   out_7507874769768704135[0] = 0;
   out_7507874769768704135[1] = 0;
   out_7507874769768704135[2] = 0;
   out_7507874769768704135[3] = 0;
   out_7507874769768704135[4] = 1;
   out_7507874769768704135[5] = 0;
   out_7507874769768704135[6] = 0;
   out_7507874769768704135[7] = 0;
   out_7507874769768704135[8] = 0;
   out_7507874769768704135[9] = 0;
   out_7507874769768704135[10] = 0;
   out_7507874769768704135[11] = 0;
   out_7507874769768704135[12] = 0;
   out_7507874769768704135[13] = 1;
   out_7507874769768704135[14] = 0;
   out_7507874769768704135[15] = 0;
}
void h_30(double *state, double *unused, double *out_7690683611848377990) {
   out_7690683611848377990[0] = state[4];
}
void H_30(double *state, double *unused, double *out_3368413205376219277) {
   out_3368413205376219277[0] = 0;
   out_3368413205376219277[1] = 0;
   out_3368413205376219277[2] = 0;
   out_3368413205376219277[3] = 0;
   out_3368413205376219277[4] = 1;
   out_3368413205376219277[5] = 0;
   out_3368413205376219277[6] = 0;
   out_3368413205376219277[7] = 0;
}
void h_26(double *state, double *unused, double *out_2689010344187532345) {
   out_2689010344187532345[0] = state[7];
}
void H_26(double *state, double *unused, double *out_6496665213197378340) {
   out_6496665213197378340[0] = 0;
   out_6496665213197378340[1] = 0;
   out_6496665213197378340[2] = 0;
   out_6496665213197378340[3] = 0;
   out_6496665213197378340[4] = 0;
   out_6496665213197378340[5] = 0;
   out_6496665213197378340[6] = 0;
   out_6496665213197378340[7] = 1;
}
void h_27(double *state, double *unused, double *out_4239011550513947736) {
   out_4239011550513947736[0] = state[3];
}
void H_27(double *state, double *unused, double *out_4655995193212844589) {
   out_4655995193212844589[0] = 0;
   out_4655995193212844589[1] = 0;
   out_4655995193212844589[2] = 0;
   out_4655995193212844589[3] = 1;
   out_4655995193212844589[4] = 0;
   out_4655995193212844589[5] = 0;
   out_4655995193212844589[6] = 0;
   out_4655995193212844589[7] = 0;
}
void h_29(double *state, double *unused, double *out_6747173758946715555) {
   out_6747173758946715555[0] = state[1];
}
void H_29(double *state, double *unused, double *out_7192211974505051230) {
   out_7192211974505051230[0] = 0;
   out_7192211974505051230[1] = 1;
   out_7192211974505051230[2] = 0;
   out_7192211974505051230[3] = 0;
   out_7192211974505051230[4] = 0;
   out_7192211974505051230[5] = 0;
   out_7192211974505051230[6] = 0;
   out_7192211974505051230[7] = 0;
}
void h_28(double *state, double *unused, double *out_1390654631818114798) {
   out_1390654631818114798[0] = state[5];
   out_1390654631818114798[1] = state[6];
}
void H_28(double *state, double *unused, double *out_3921036825704352615) {
   out_3921036825704352615[0] = 0;
   out_3921036825704352615[1] = 0;
   out_3921036825704352615[2] = 0;
   out_3921036825704352615[3] = 0;
   out_3921036825704352615[4] = 0;
   out_3921036825704352615[5] = 1;
   out_3921036825704352615[6] = 0;
   out_3921036825704352615[7] = 0;
   out_3921036825704352615[8] = 0;
   out_3921036825704352615[9] = 0;
   out_3921036825704352615[10] = 0;
   out_3921036825704352615[11] = 0;
   out_3921036825704352615[12] = 0;
   out_3921036825704352615[13] = 0;
   out_3921036825704352615[14] = 1;
   out_3921036825704352615[15] = 0;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_312213572642374608) {
  err_fun(nom_x, delta_x, out_312213572642374608);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_486699098684977001) {
  inv_err_fun(nom_x, true_x, out_486699098684977001);
}
void car_H_mod_fun(double *state, double *out_3333653717063534221) {
  H_mod_fun(state, out_3333653717063534221);
}
void car_f_fun(double *state, double dt, double *out_1478154795453816204) {
  f_fun(state,  dt, out_1478154795453816204);
}
void car_F_fun(double *state, double dt, double *out_2267919153160654385) {
  F_fun(state,  dt, out_2267919153160654385);
}
void car_h_25(double *state, double *unused, double *out_2116182985142631742) {
  h_25(state, unused, out_2116182985142631742);
}
void car_H_25(double *state, double *unused, double *out_6245485705572964003) {
  H_25(state, unused, out_6245485705572964003);
}
void car_h_24(double *state, double *unused, double *out_5094452919565144864) {
  h_24(state, unused, out_5094452919565144864);
}
void car_H_24(double *state, double *unused, double *out_7507874769768704135) {
  H_24(state, unused, out_7507874769768704135);
}
void car_h_30(double *state, double *unused, double *out_7690683611848377990) {
  h_30(state, unused, out_7690683611848377990);
}
void car_H_30(double *state, double *unused, double *out_3368413205376219277) {
  H_30(state, unused, out_3368413205376219277);
}
void car_h_26(double *state, double *unused, double *out_2689010344187532345) {
  h_26(state, unused, out_2689010344187532345);
}
void car_H_26(double *state, double *unused, double *out_6496665213197378340) {
  H_26(state, unused, out_6496665213197378340);
}
void car_h_27(double *state, double *unused, double *out_4239011550513947736) {
  h_27(state, unused, out_4239011550513947736);
}
void car_H_27(double *state, double *unused, double *out_4655995193212844589) {
  H_27(state, unused, out_4655995193212844589);
}
void car_h_29(double *state, double *unused, double *out_6747173758946715555) {
  h_29(state, unused, out_6747173758946715555);
}
void car_H_29(double *state, double *unused, double *out_7192211974505051230) {
  H_29(state, unused, out_7192211974505051230);
}
void car_h_28(double *state, double *unused, double *out_1390654631818114798) {
  h_28(state, unused, out_1390654631818114798);
}
void car_H_28(double *state, double *unused, double *out_3921036825704352615) {
  H_28(state, unused, out_3921036825704352615);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);
