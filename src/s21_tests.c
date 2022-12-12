#include <check.h>
#include <math.h>
#include <stdlib.h>

#include "s21_math.h"

START_TEST(abs_test) {
  ck_assert_int_eq(abs(-15), s21_abs(-15));
  ck_assert_int_eq(abs(15), s21_abs(15));
  ck_assert_int_eq(abs(-100), s21_abs(-100));
  ck_assert_int_eq(abs(45743785), s21_abs(45743785));
  ck_assert_int_eq(abs(0), s21_abs(0));
}
END_TEST

START_TEST(fabs_test) {
  ck_assert_int_eq(fabs(-15.468421), s21_fabs(-15.468421));
  ck_assert_int_eq(fabs(15.468421), s21_fabs(15.468421));
  ck_assert_int_eq(fabs(-1000.23453), s21_fabs(-1000.23453));
  ck_assert_int_eq(fabs(45743785.297802), s21_fabs(45743785.297802));
  ck_assert_int_eq(fabs(0.0), s21_fabs(0.0));
}
END_TEST

START_TEST(fmod_test) {
  ck_assert_ldouble_eq(fmod(1000.23453, 3.2), s21_fmod(1000.23453, 3.2));
  ck_assert_ldouble_eq(fmod(5743.297802, 1.4), s21_fmod(5743.297802, 1.4));
  ck_assert_ldouble_eq(fmod(0, 5.3), s21_fmod(0, 5.3));
  ck_assert_ldouble_eq(fmod(15.468421, 8.8), s21_fmod(15.468421, 8.8));
  ck_assert_ldouble_eq(fmod(13.21313, 4.67), s21_fmod(13.21313, 4.67));
  ck_assert_ldouble_eq(fmod(14.1, 14.1), s21_fmod(14.1, 14.1));
}
END_TEST

START_TEST(pow_test) {
  ck_assert_float_eq(pow(45743785.297802, 1.4), s21_pow(45743785.297802, 1.4));
  ck_assert_float_eq(pow(0, 5.3), s21_pow(0, 5.3));
  ck_assert_float_eq(pow(15.468421, 8.8), s21_pow(15.468421, 8.8));
  ck_assert_int_eq(pow(-1000.23453, 3.2), s21_pow(-1000.23453, 3.2));
  ck_assert_float_eq(pow(562, -4.0), s21_pow(562, -4.0));
  ck_assert_int_eq(pow(-15.468421, 8.8), s21_pow(-15.468421, 8.8));
}
END_TEST

START_TEST(ceil_test) {
  ck_assert_ldouble_eq(ceil(1000.0000001), s21_ceil(1000.0000001));
  ck_assert_ldouble_eq(ceil(15.468421), s21_ceil(15.468421));
  ck_assert_ldouble_eq(ceil(3000), s21_ceil(3000));
  ck_assert_ldouble_eq(ceil(-15.468421), s21_ceil(-15.468421));
  ck_assert_ldouble_eq(ceil(0), s21_ceil(0));
  ck_assert_ldouble_eq(ceil(45743785.297802), s21_ceil(45743785.297802));
}
END_TEST

START_TEST(floor_test) {
  ck_assert_ldouble_eq(floor(1000.0000001), s21_floor(1000.0000001));
  ck_assert_ldouble_eq(floor(-15.468), s21_floor(-15.468));
  ck_assert_ldouble_eq(floor(-0.231), s21_floor(-0.231));
  ck_assert_ldouble_eq(floor(3000), s21_floor(3000));
  ck_assert_ldouble_eq(floor(4.67), s21_floor(4.67));
  ck_assert_ldouble_eq(floor(0.0), s21_floor(0.0));
}
END_TEST

START_TEST(exp_test) {
  ck_assert_int_eq(exp(s21_nan), s21_exp(s21_nan));
  ck_assert_float_eq(exp(s21_DBL_MAX), s21_exp(s21_DBL_MAX));
  ck_assert_float_eq(exp(1), s21_exp(1));
  ck_assert_float_eq(exp(-1), s21_exp(-1));
  ck_assert_double_eq(exp(100), s21_exp(100));
  ck_assert_float_eq(exp(-1), s21_exp(-1));
  ck_assert_float_eq(exp(0), s21_exp(0));
  ck_assert_float_eq(exp(15.468421), s21_exp(15.468421));
  ck_assert_float_eq(exp(-15.468421), s21_exp(-15.468421));
  ck_assert_float_eq(exp(3000), s21_exp(3000));
  ck_assert_float_eq(exp(0.231), s21_exp(0.231));
  ck_assert_double_eq(exp(1000.0000001), s21_exp(1000.0000001));
}
END_TEST

START_TEST(sin_test) {
  ck_assert_float_eq(sin(0.0), s21_sin(0.0));
  ck_assert_float_eq(sin(1.0), s21_sin(1.0));
  ck_assert_float_eq(sin(-1.0), s21_sin(-1.0));
  ck_assert_float_eq(sin(-15.4684), s21_sin(-15.4684));
  ck_assert_float_eq(sin(15.4684), s21_sin(15.4684));
  ck_assert_float_eq(sin(0.342), s21_sin(0.342));
  ck_assert_float_eq(sin(-0.98), s21_sin(-0.98));
}
END_TEST

START_TEST(cos_test) {
  ck_assert_float_eq(sin(0.0), s21_sin(0.0));
  ck_assert_float_eq(sin(1.0), s21_sin(1.0));
  ck_assert_float_eq(sin(-1.0), s21_sin(-1.0));
  ck_assert_float_eq(sin(-15.4684), s21_sin(-15.4684));
  ck_assert_float_eq(sin(15.4684), s21_sin(15.4684));
  ck_assert_float_eq(sin(0.342), s21_sin(0.342));
  ck_assert_float_eq(sin(-0.98), s21_sin(-0.98));
}
END_TEST

START_TEST(tan_test) {
  ck_assert_float_eq(sin(0.0), s21_sin(0.0));
  ck_assert_float_eq(sin(1.0), s21_sin(1.0));
  ck_assert_float_eq(sin(-1.0), s21_sin(-1.0));
  ck_assert_float_eq(sin(-15.4684), s21_sin(-15.4684));
  ck_assert_float_eq(sin(15.4684), s21_sin(15.4684));
  ck_assert_float_eq(sin(0.342), s21_sin(0.342));
  ck_assert_float_eq(sin(-0.98), s21_sin(-0.98));
  ck_assert_float_eq(sin(155.0), s21_sin(155.0));
}
END_TEST

START_TEST(sqrt_test) {
  ck_assert_float_eq(sqrt(0.0), s21_sqrt(0.0));
  ck_assert_float_eq(sqrt(1.0), s21_sqrt(1.0));
  ck_assert_float_eq(sin(-15.4684), s21_sin(-15.4684));
  ck_assert_float_eq(sin(15.4684), s21_sin(15.4684));
  ck_assert_float_eq(sin(0.342), s21_sin(0.342));
  ck_assert_float_eq(sin(-0.98), s21_sin(-0.98));
  ck_assert_float_eq(sin(155.0), s21_sin(155.0));
}

START_TEST(log_test) {
  ck_assert_float_eq(log(0.0), s21_log(0.0));
  ck_assert_int_eq(log(1.0), s21_log(1.0));
  ck_assert_int_eq(log(-15.4684), s21_log(-15.4684));
  ck_assert_float_eq(log(15.4684), s21_log(15.4684));
  ck_assert_float_eq(log(2.245), s21_log(2.245));
  ck_assert_float_eq(log(43785.297), s21_log(43785.297));
  // ck_assert_float_eq(log(-0.98), s21_log(-0.98)); ?
}
END_TEST

START_TEST(asin_test) {
  ck_assert_float_eq(asin(0), s21_asin(0));
  ck_assert_float_eq(asin(1), s21_asin(1));
  ck_assert_float_eq(asin(-1), s21_asin(-1));
  ck_assert_int_eq(asin(-15.4684), s21_asin(-15.4684));
  ck_assert_int_eq(asin(15.4684), s21_asin(15.4684));
  ck_assert_float_eq(asin(0.98), s21_asin(0.98));
  ck_assert_float_eq(asin(0.5), s21_asin(0.5));
  ck_assert_float_eq(asin(-0.5), s21_asin(-0.5));
  ck_assert_float_eq(asin(s21_EPS), s21_asin(s21_EPS));
}
END_TEST

START_TEST(acos_test) {
  ck_assert_float_eq(acos(0), s21_acos(0));
  ck_assert_float_eq(acos(1), s21_acos(1));
  ck_assert_float_eq(acos(-1), s21_acos(-1));
  ck_assert_int_eq(acos(-15.4684), s21_acos(-15.4684));
  ck_assert_int_eq(acos(15.4684), s21_acos(15.4684));
  ck_assert_float_eq(acos(0.98), s21_acos(0.98));
  ck_assert_float_eq(acos(0.5), s21_acos(0.5));
  ck_assert_float_eq(acos(-0.5), s21_acos(-0.5));
}
END_TEST

START_TEST(atan_test) {
  ck_assert_float_eq(atan(0), s21_atan(0));
  ck_assert_float_eq(atan(1), s21_atan(1));
  ck_assert_float_eq(atan(-1), s21_atan(-1));
  ck_assert_int_eq(atan(-15.4684), s21_atan(-15.4684));
  ck_assert_float_eq(atan(15.4684), s21_atan(15.4684));
  ck_assert_float_eq(atan(3.14), s21_atan(3.14));
  ck_assert_float_eq(atan(4.67), s21_atan(4.67));
  ck_assert_float_eq(atan(2.2343), s21_atan(2.2343));
  ck_assert_float_eq(atan(0.98), s21_atan(0.98));
  // ck_assert_float_eq(atan(-0.98), s21_atan(-0.98)); ?
  ck_assert_float_eq(atan(43785.297), s21_atan(43785.297));
}
END_TEST

Suite *s21_math_tests_create() {
  Suite *s21_math = suite_create("s21_math");
  TCase *s21_math_tests = tcase_create("S21_MATH");
  tcase_add_test(s21_math_tests, ceil_test);
  tcase_add_test(s21_math_tests, floor_test);
  tcase_add_test(s21_math_tests, exp_test);
  tcase_add_test(s21_math_tests, sin_test);
  tcase_add_test(s21_math_tests, cos_test);
  tcase_add_test(s21_math_tests, tan_test);
  tcase_add_test(s21_math_tests, abs_test);
  tcase_add_test(s21_math_tests, fabs_test);
  tcase_add_test(s21_math_tests, fmod_test);
  tcase_add_test(s21_math_tests, pow_test);
  tcase_add_test(s21_math_tests, sqrt_test);
  tcase_add_test(s21_math_tests, log_test);
  tcase_add_test(s21_math_tests, asin_test);
  tcase_add_test(s21_math_tests, acos_test);
  tcase_add_test(s21_math_tests, atan_test);
  suite_add_tcase(s21_math, s21_math_tests);
  return s21_math;
}

int main() {
  Suite *s21_math = s21_math_tests_create();
  SRunner *s21_math_runner = srunner_create(s21_math);
  int number_failed;
  srunner_run_all(s21_math_runner, CK_NORMAL);
  number_failed = srunner_ntests_failed(s21_math_runner);
  srunner_free(s21_math_runner);
  return number_failed == 0 ? 0 : 1;
}
