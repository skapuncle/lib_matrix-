CC=g++
LDFLAGS=-lgtest
CFLAGS=-Wall -Werror -Wextra -std=c++17
SRC=s21_matrix_oop.cc
TEST=tests.cc
TESTO=$(TEST:%.cc=%.o)
TEST_OBJ=$(SRC:%.cc=%.o)
LIB=s21_matrix_oop.a
LIBO=$(LIB:%.a=%.o)

all: test gcov_report

test: clean $(TEST_OBJ) $(TESTO)
	@${CC} ${TEST_OBJ} ${TESTO} ${CFLAGS} ${LDFLAGS} -o test --coverage
	@./test

%.o: %.cc
	@${CC} ${CFLAGS} -c $^

s21_matrix_oop.a: s21_matrix_oop.h $(LIBO)
	@ar rcs $@ $^
	@ranlib $@

$(TEST_OBJ): $(SRC)
	@$(CC) $(CFLAGS) -c $^ --coverage

clean:
	@/bin/rm -rf *.o *.a test *.gcno *gcda report *.info test *.out

gcov_report: test
	@lcov -t "test" -o test.info --no-external -c -d .
	@genhtml -o report test.info
	@open report/index.html

style:
	clang-format --style=google -i *.cc *.h

cppcheck:
	cppcheck --enable=all --suppress=missingIncludeSystem .

valgrind:
	valgrind --trace-children=yes --track-fds=yes --track-origins=yes --leak-check=full --show-leak-kinds=all ./test
