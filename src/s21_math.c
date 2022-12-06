#include <stdio.h>
#include <math.h>
#include "s21_math.h"

int s21_abs(int x) {
    return x < 0 ? -x : x;
}

long double s21_fabs(double x) {
    return x < 0 ? -x : x;
}

int s21_factorial (int n) {
  return (n < 2) ? 1 : n * s21_factorial (n - 1);
}

long double s21_sin(double x) {
    for (; x < -2 * s21_PI || 2 * s21_PI < x;) {
        if (x > 2 * s21_PI) {
            x -= 2 * s21_PI;
        } else {
            x += 2 * s21_PI;
        }
    }
    double result = x;
    double elem = x;
    long double ch = x;
    long double zn = 1.0;
    for (int i = 1; i <= 10000; i++) {
            ch = x * x;
            zn = (2.0 * i) * ((2.0 * i) + 1);
            elem = elem * ch/zn;
        if (i % 2 == 0)
            result += elem;
        else
            result -= elem;
    }
    return result;
}

long double s21_cos(double x) {
    for (; x < -2 * s21_PI || 2 * s21_PI < x;) {
        if (x > 2 * s21_PI) {
            x -= 2 * s21_PI;
        } else {
            x += 2 * s21_PI;
        }
    }
    double result = 1;
    double elem = 1;
    long double ch = 1;
    long double zn = 1.0;
    for (int i = 1; i <= 10000; i++) {
            ch = x * x;
            zn = (2.0 * i) * ((2.0 * i) - 1);
            elem = elem * ch/zn;
        if (i % 2 == 0)
            result += elem;
        else
            result -= elem;
    }
    return result;
}

long double s21_tan(double x) {
    double tang = s21_sin(x) / s21_cos(x);
    return tang;
}

long double s21_atan(double x) {
    long double res = x;
    long double ch = x;
    long double zn = 1;
    long double elem = x;
    int flag = 0;
    int sign = 1;
    if (x < 0) {
        sign = -1;
        x = fabs(x);
    }
    if (x > 1 || x < -1) {
        flag = 1;
        ch = 1 / x;
        res = 1 / x;
        x = 1 / x;
    }
    for (int n = 1; n < 500; n++) {
        ch = ch * x * x;
        zn = n * 2 + 1;
        elem = ch / zn;
        if (n % 2 == 0)
            res += elem;
        else
            res -= elem;
    }
    if (flag) {
        res = s21_PI / 2 - res;
    }
    if (x == 1) {
        res = s21_PI / 4;
    }
    return sign * res;
}

long double s21_sqrt(double x) {
    long double s1 = 0, s2 = 4;
    while (s21_fabs(s2 - s1) > 1e-10) {
        if (x < 0) {
            break;
        }
        s1 = s2;
        s2 = (s1 + x / s1) / 2;
    }
    return s2;
}

long double s21_asin(double x) {
    long double result = 0;
    int sign = 1;
    if (x < 0) {
        sign = -1;
        x = fabs(x);
    }
    if (x == 1) {
        result = s21_PI / 2;
    } else if (x == -1) {
        result = -s21_PI / 2;
    } else if (x < 1 && x > -1) {
        result = s21_atan(x / s21_sqrt(1 - x * x));
    }
    return sign * result;
}

long double s21_acos(double x) {
    return (x <= 1 && x >= -1) ? s21_PI / 2 - s21_asin(x) : 0.;
}

long double s21_log(double x) {
    return (x > 1) ? 1 + s21_log(x / 5) : 0;
}

long double s21_exp(double x) {
    long double res = 1, y = 1, n = 1;
    for (; s21_fabs((double)res) > s21_EPS; y+=res) {
        if (y > s21_DBL_MAX || y < s21_DBL_MIN) {
            y = (y > s21_DBL_MAX) ? s21_inf : 0.;
            break;
        }
        res *= x/n++;
    }
    return y;
}

long double s21_fmod(double x, double y) {
    long double rem = s21_fabs(x);
    for(int i = 1; rem > s21_fabs(y); i++) {
        rem = s21_fabs(x) - y * i;
    }
    return x < 0 ? -rem : rem;
}

int main() {
    printf("kon = %lf\n", fmod(-234.432, 23));
    printf("u nas = %Lf\n", s21_fmod(-234.432,23));
}

