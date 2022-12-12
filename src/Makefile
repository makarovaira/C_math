CFLAGS = -Wall -Werror -Wextra -std=c11
GCOV=-fprofile-arcs -ftest-coverage
CC = gcc
LIBS = -lcheck
RM = rm -rf
SRCS = s21_math.c
OBJS = $(SRCS:.c=.o)

all: s21_math.a test gcov_report

s21_math.a:
	$(CC) -c $(CFLAGS) $(SRCS)
	ar rcs s21_math.a $(OBJS)
	ranlib s21_math.a

test:
	$(CC) s21_tests.c -lcheck -L. s21_math.a -o s21_test
	./s21_test

gcov_report:
	$(CC) $(GCOV) s21_tests.c $(SRCS) -o s21_test -lcheck
	./s21_test
	lcov -t "test" -o test.info -c -d ./
	genhtml test.info -o report
	open report/index.html

create_test:
	checkmk clean_mode=1 s21_test_math.check > s21_test_math.c

install_brew:
	cd ~
	curl -fsSL https://rawgit.com/kube/42homebrew/master/install.sh | zsh
	brew install lcov

clean:
	$(RM) *.a *.o *.info *.gcda *.gcno *.gcov s21_test *.a

rebuild:
	$(MAKE) clean
	$(MAKE) all
