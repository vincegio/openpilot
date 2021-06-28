#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_312213572642374608);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_486699098684977001);
void car_H_mod_fun(double *state, double *out_3333653717063534221);
void car_f_fun(double *state, double dt, double *out_1478154795453816204);
void car_F_fun(double *state, double dt, double *out_2267919153160654385);
void car_h_25(double *state, double *unused, double *out_2116182985142631742);
void car_H_25(double *state, double *unused, double *out_6245485705572964003);
void car_h_24(double *state, double *unused, double *out_5094452919565144864);
void car_H_24(double *state, double *unused, double *out_7507874769768704135);
void car_h_30(double *state, double *unused, double *out_7690683611848377990);
void car_H_30(double *state, double *unused, double *out_3368413205376219277);
void car_h_26(double *state, double *unused, double *out_2689010344187532345);
void car_H_26(double *state, double *unused, double *out_6496665213197378340);
void car_h_27(double *state, double *unused, double *out_4239011550513947736);
void car_H_27(double *state, double *unused, double *out_4655995193212844589);
void car_h_29(double *state, double *unused, double *out_6747173758946715555);
void car_H_29(double *state, double *unused, double *out_7192211974505051230);
void car_h_28(double *state, double *unused, double *out_1390654631818114798);
void car_H_28(double *state, double *unused, double *out_3921036825704352615);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}