#include "s21_math.h"

int s21_abs(int x) { return x < 0 ? -x : x; }

long double s21_fabs(double x) { return x < 0 ? -x : x; }

long double s21_sin(double x) {
  double result, elem;
  if (s21_fabs(x) == s21_inf) {
    result = -s21_nan;
  } else {
    x = (x > 2 * s21_PI || x < -2 * s21_PI) ? s21_fmod(x, 2 * s21_PI) : x;
    result = elem = x;
    for (int i = 1; s21_fabs(elem) > s21_EPS; i++) {
      elem *= x * x / ((2.0 * i) * ((2.0 * i) + 1));
      result += i % 2 ? -elem : elem;
    }
  }
  return result;
}

long double s21_cos(double x) {
  return s21_fabs(x) == s21_inf ? -s21_nan : s21_sin(s21_PI / 2. - x);
}

long double s21_tan(double x) { return s21_sin(x) / s21_cos(x); }

long double s21_atan(double x) {
  long double res, ch, elem;
  int flag_span = 0, sign = (x < 0) ? 0 : 1;
  x = s21_fabs(x);
  res = ch = elem = x;
  if (x > 1) {
    flag_span = 1;
    ch = res = x = 1 / x;
  }
  if (x == 1) {
    res = s21_PI / 4;
  } else {
    for (int n = 1; s21_fabs(elem) > s21_EPS; n++) {
      ch *= x * x;
      elem = ch / (n * 2 + 1);
      res += n % 2 ? -elem : elem;
    }
  }
  res = (flag_span) ? s21_PI / 2 - res : res;
  return (sign) ? res : -res;
}

long double s21_sqrt(double x) {
  long double rez = 4, y = 0;
  if (x == s21_inf) {
    rez = s21_inf;
  } else if (x == -s21_inf) {
    rez = s21_nan;
  } else {
    while (s21_fabs(rez - y) > s21_EPS) {
      if (x < 0) {
        rez = s21_nan;
        break;
      }
      y = rez;
      rez = (y + x / y) / 2;
    }
  }
  return x == 0 ? 0 : rez;
}

long double s21_asin(double x) {
  long double result;
  int sign = (x < 0) ? 0 : 1;
  x = s21_fabs(x);
  if (x == 1) {
    result = s21_PI / 2;
  } else if (x < 1 && x >= 0) {
    result = s21_atan(x / s21_sqrt(1 - x * x));
  } else {
    result = s21_nan;
    x == s21_inf ? sign = 0 : sign;
  }
  return (sign) ? result : -result;
}

long double s21_acos(double x) {
  return (x <= 1 && x >= -1)      ? s21_PI / 2 - s21_asin(x)
         : s21_fabs(x) == s21_inf ? -(s21_nan)
                                  : s21_nan;
}

long double s21_exp(double x) {
  long double res = 1, y = 1, n = 1;
  int neg = (x < 0) ? 1 : 0;
  x = s21_fabs(x);
  while (s21_fabs((double)res) > s21_EPS) {
    if (y > s21_DBL_MAX) {
      y = s21_inf;
      break;
    }
    res *= x / n++;
    y += res;
  }
  y = (neg) ? y > s21_DBL_MAX ? 0 : 1. / y : y;
  return y;
}

long double s21_log(double x) {
  long double res;
  if (x == s21_inf) {
    res = s21_inf;
  } else if (x == -s21_inf) {
    res = s21_nan;
  } else {
    unsigned a = 0, d = 1;
    long double b = 0, c = (x < 1) ? 1 / x : x, e = b, f;
    if (x > 0) {
      for (; (c /= s21_EXP) > 1; ++a)
        ;
      c = 2 * (1 / (c * s21_EXP - 1)) + 1;
      f = c * c;
      for (c /= 2; b += 1 / (d * c), b - e; d += 2) {
        c *= f;
        e = b;
      }
    } else
      b = (x == 0) / 0.;
    res = a + b;
  }
  return s21_fabs(x) == s21_inf ? res : x < 1 ? -res : res;
}

long double s21_fmod(double x, double y) {
  long double rem;
  if (x == s21_inf || y == 0) {
    rem = s21_nan;
  } else if (y == s21_inf || y == -s21_inf) {
    rem = (long long)x;
  } else {
    long double dx = x;
    long double dy = y;
    long long div = dx / dy;
    dx -= div * dy;
    rem = dx;
  }
  return rem;
}

long double s21_pow(double base, double exp) {
  long double power;
  if ((long int)exp == exp) {
    if (exp > 0) {
      power = base;
      for (long int i = 0; i < (long int)exp - 1; i++) {
        power *= base;
      }
    } else if (exp == 0) {
      power = 1;
    } else {
      power = 1 / base;
      for (long int i = 0; i < (long int)exp * (-1) - 1; i++) {
        power /= base;
      }
    }
  } else if (base < 0) {
    if (s21_fabs(exp) == s21_inf) {
      if (s21_fabs(base) < 1) {
        power = 0;
      } else if (s21_fabs(base) == 1) {
        power = 1;
      } else {
        power = (exp == -s21_inf) ? 0 : s21_inf;
      }
    } else {
      power = -s21_nan;
    }
  } else if (base == 0) {
    power = (exp == 0) ? 1 : 0;
  } else if (base == 1) {
    power = 1;
  } else {
    power = s21_exp(exp * (double)s21_log(base));
  }
  return power;
}

long double s21_floor(double x) {
  long double res = x;
  if (s21_fabs(x) == s21_inf) {
    res = x;
  } else if (x == s21_nan) {
    res = s21_nan;
  } else {
    if (res >= 0) {
      res -= s21_fmod(x, 1);
    } else if (s21_fabs(s21_fmod(x, 1)) > s21_EPS) {
      res -= 1;
      res -= s21_fmod(x, 1);
    }
  }
  return res;
}

long double s21_ceil(double x) {
  long double res = x;
  if (s21_fabs(x) == s21_inf) {
    res = x;
  } else if (x == s21_nan) {
    res = s21_nan;
  } else if (s21_fabs(s21_fmod(x, 1)) > 0) {
    if (x > 0.0) {
      res = (long int)(x + 1.0);
    } else {
      res = (long int)x;
    }
  }
  return res;
}
