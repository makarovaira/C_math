#ifndef SRC_S21_MATH_H_
#define SRC_S21_MATH_H_

#define s21_PI 3.14159265358979323846264338327950288
#define s21_EPS 1e-17
#define s21_EXP 2.7182818284590452353602874713526624
#define s21_E 2.71828182845904523536028747135266249
#define s21_DBL_MAX 1.7976931348623157e308L
#define s21_DBL_MIN 2.2250738585072014e-308L
#define s21_inf 1.0 / 0.0
#define s21_nan 0.0 / 0.0

int s21_abs(int x);
long double s21_fabs(double x);
long double s21_sin(double x);
long double s21_cos(double x);
long double s21_tan(double x);
long double s21_atan(double x);
long double s21_acos(double x);
long double s21_sqrt(double x);
long double s21_asin(double x);
long double s21_exp(double x);
long double s21_log(double x);
long double s21_fmod(double x, double y);
long double s21_pow(double base, double exp);
long double s21_floor(double x);
long double s21_ceil(double x);

#endif  // SRC_S21_MATH_H_
