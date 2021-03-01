CC = g++
CF =

CPPFILES = $(shell find src -name '*.cpp')
HPPFILES = $(shell find headers -name '*.hpp')
OFILES = $(addprefix obj/, $(CPPFILES:src/%.c=%.o))

all: obj libSparseMatrix.a

obj:
	mkdir -p obj

obj/%.o: src/%.cpp $(HPPFILES)
	@$(CC) $(CF) -c $< -o $@ -I headers

libSparseMatrix.a: $(OFILES)
	@ar rc $@ $(OFILES)
	@ranlib $@

clean:
	rm -rf obj

fclean: clean
	rm -f SparseMatrix.a

re: fclean all

.PHONY: all clean fclean re
