#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_3(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_19(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_7100013140337143968);
void live_err_fun(double *nom_x, double *delta_x, double *out_6609159293288784210);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_8190965837586844167);
void live_H_mod_fun(double *state, double *out_9091563194564624183);
void live_f_fun(double *state, double dt, double *out_2084886017125782583);
void live_F_fun(double *state, double dt, double *out_7865388672361347508);
void live_h_3(double *state, double *unused, double *out_5721050071790325762);
void live_H_3(double *state, double *unused, double *out_4896709827861440587);
void live_h_4(double *state, double *unused, double *out_6958619119480655761);
void live_H_4(double *state, double *unused, double *out_3961471467626008537);
void live_h_9(double *state, double *unused, double *out_5982385983887379598);
void live_H_9(double *state, double *unused, double *out_6283977011496864507);
void live_h_10(double *state, double *unused, double *out_6316084181607648543);
void live_H_10(double *state, double *unused, double *out_1030569411575115501);
void live_h_12(double *state, double *unused, double *out_3629879685827352034);
void live_H_12(double *state, double *unused, double *out_6068026979365218106);
void live_h_31(double *state, double *unused, double *out_1743096937849468433);
void live_H_31(double *state, double *unused, double *out_3392937204221798285);
void live_h_32(double *state, double *unused, double *out_4966958122310810386);
void live_H_32(double *state, double *unused, double *out_2907428540038602901);
void live_h_13(double *state, double *unused, double *out_2080022393384731655);
void live_H_13(double *state, double *unused, double *out_8127275896522900919);
void live_h_14(double *state, double *unused, double *out_5982385983887379598);
void live_H_14(double *state, double *unused, double *out_6283977011496864507);
void live_h_19(double *state, double *unused, double *out_2180750534064810423);
void live_H_19(double *state, double *unused, double *out_684057922475399103);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}