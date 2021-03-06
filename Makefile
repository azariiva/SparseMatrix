
CC = g++
CF = -std=c++11 -g -ftrapv -fsanitize=undefined -Wall -Werror -Wformat-security -Wignored-qualifiers -Winit-self -Wswitch-default -Wfloat-equal -Wshadow -Wpointer-arith -Wtype-limits -Wempty-body -Wlogical-op -Wmissing-field-initializers -Wcast-qual -Wwrite-strings -lm -L.

CPPFILES = $(shell find src -name '*.cpp')
HPPFILES = $(shell find headers -name '*.hpp')
OFILES = $(addprefix obj/, $(CPPFILES:src/%.cpp=%.o))

CPPTESTS =  $(shell find tests -name '*.cpp')
OUTTESTS = $(CPPTESTS:%.cpp=%.out)

all: obj libSparseMatrix.a

obj:
	mkdir -p obj

obj/%.o: src/%.cpp $(HPPFILES)
	@$(CC) $(CF) -I headers -c $< -o $@
	@echo "$@ compiled"

libSparseMatrix.a: $(OFILES)
	@ar rc $@ $(OFILES)
	@ranlib $@
	@echo "$@ compiled"

clean:
	rm -rf obj

fclean: clean
	rm -f libSparseMatrix.a
	rm -f $(OUTTESTS)

re: fclean all

tests: $(OUTTESTS)

tests/%.out: tests/%.cpp all
	$(CC) $(CF) -I headers -L. $< -lSparseMatrix -o $@

.PHONY: all clean fclean re tests
