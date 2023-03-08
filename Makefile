include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
#CC=gcc
#FC = ifort
CLINKER := mpicc -fPIC #-wd1572 -O3 -axAVX,SSE4.2 -fno-alias -no-prec-div -no-prec-sqrt -ip
MAKE    := /usr/bin/make
MKDIR_P := /bin/mkdir -p
OBJ_DIR := obj
LIB_DIR := libs
BIN_DIR := bin
EXE     := ${BIN_DIR}/ex1
SRC_DIR := src
SOURCES := $(shell find . -type f -name '*.c')
HEADERS := $(shell find . -type f -name '*.h')
OBJ     := $(addsuffix .o,$(basename $(SOURCES)))
DEP     := $(OBJ:.o=.d)
# Preprocessor flags here
CPPFLAGS:=  -MMD -MP -I. 
# Compiler flags here
CFLAGS  := -std=c99  ${SLEPC_INCLUDE} ${PETSC_CC_INCLUDES} 
LDFLAGS := ${SLEPC_INCLUDE} ${PETSC_CC_INCLUDES} 
LDLIBS  := ${LIB_DIR}/irpf90.a ${SLEPC_EPS_LIB}
 
.PHONY: ex1 all debug run clean test

ex1: ${BIN_DIR}/ex1

${OBJ_DIR}:
	${MKDIR_P} ${OBJ_DIR}

${LIB_DIR}:
	${MKDIR_P} ${LIB_DIR}

${BIN_DIR}:
	${MKDIR_P} ${BIN_DIR}

directories: ${OBJ_DIR} ${LIB_DIR} ${BIN_DIR}

${LIB_DIR}/irpf90.a: directories
	cd ${SRC_DIR} && irpf90 init && $(MAKE) && cp IRPF90_temp/irpf90.a ../${LIB_DIR}


all: CPPFLAGS += -DNDEBUG
all: $(EXE)

$(EXE): $(OBJ)
	@echo "Compilation complete!"
	$(CLINKER) $(LDFLAGS) $^ $(LDLIBS) -o $@ 
	@echo "Linking complete!"

-include $(DEP)

test:
	@echo ${SOURCES}
	@echo ${HEADERS}
	#@echo ${OBJ}
	#@echo ${LDFLAGS}

