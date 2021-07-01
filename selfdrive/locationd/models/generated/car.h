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
void car_err_fun(double *nom_x, double *delta_x, double *out_6174462923116197084);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_1639972823104049711);
void car_H_mod_fun(double *state, double *out_3507200293686001885);
void car_f_fun(double *state, double dt, double *out_4902583061038932869);
void car_F_fun(double *state, double dt, double *out_7818687564155864474);
void car_h_25(double *state, double *unused, double *out_7190281350844188701);
void car_H_25(double *state, double *unused, double *out_4506058534437431261);
void car_h_24(double *state, double *unused, double *out_3285867608003353900);
void car_H_24(double *state, double *unused, double *out_6898639677699028317);
void car_h_30(double *state, double *unused, double *out_494511361579178926);
void car_H_30(double *state, double *unused, double *out_4326786628322937075);
void car_h_26(double *state, double *unused, double *out_673232439496962811);
void car_H_26(double *state, double *unused, double *out_1198534620501778012);
void car_h_27(double *state, double *unused, double *out_4585876030035491050);
void car_H_27(double *state, double *unused, double *out_3039204640486311763);
void car_h_29(double *state, double *unused, double *out_3307063464849032309);
void car_H_29(double *state, double *unused, double *out_7549017147828961947);
void car_h_28(double *state, double *unused, double *out_4201569591468013248);
void car_H_28(double *state, double *unused, double *out_4530271917774028433);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}