SRCDIR = ./src
OBJDIR = ./obj
BINDIR = ./bin
EXEC = linalg

SRC = $(wildcard $(SRCDIR)/*.cpp)
OBJ = $(addprefix $(OBJDIR)/,$(notdir $(SRC:.cpp=.o)))

CC          = g++-4.7
CFLAGS      = -c -Wall -std=c++11 -Os
LDFLAGS     = -lm -lstdc++

X86_64_FLAGS = -mmmx -msse -msse2 -m3dnow -m64 -malign-double -mieee-fp -O2
ARM_FLAGS    = -mfloat-abi=hard -mfpu=neon -march=armv7-a -mtune=cortex-a8 -ffast-math -O3

ARCH = $(shell lscpu | grep Architecture | awk '{print $$2}')

ifeq (${ARCH}, x86_64)
	CFLAGS += $(X86_64_FLAGS)
endif

ifeq (${ARCH}, armv7l)
   CFLAGS += $(ARM_FLAGS)
endif


.PHONY: all
all: compile link

.PHONY: compile
compile: $(OBJ)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) $< -o $@

.PHONY: link
link: compile
	$(CC) -o $(EXEC) $(LDFLAGS) $(OBJ)

.PHONY: test
test:
	./test.sh

.PHNOY : clean
clean:
	find -type f -name "$(EXEC)" -delete
	find -type f -name "*.o" -delete
	find -type f -name "*~" -delete

