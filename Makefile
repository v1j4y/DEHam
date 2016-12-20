include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
#CC=gcc
#FC = ifort
MAKE = /usr/bin/make
MKDIR_P = /bin/mkdir -p
OBJ_DIR := obj
LIB_DIR := libs
BIN_DIR := bin
SRC_DIR := src

.PHONY: ex1

ex1: ${BIN_DIR}/ex1

${OBJ_DIR}:
	${MKDIR_P} ${OBJ_DIR}

${LIB_DIR}:
	${MKDIR_P} ${LIB_DIR}

${BIN_DIR}:
	${MKDIR_P} ${BIN_DIR}

directories: ${OBJ_DIR} ${LIB_DIR} ${BIN_DIR}

${LIB_DIR}/irpf90.a: directories
	cd ${SRC_DIR} && irpf90 init && $(MAKE) irpf90.a && cp irpf90.a ../${LIB_DIR}

${OBJ_DIR}/read2.o: ${SRC_DIR}/read2.c directories chkopts
	${CC} ${SLEPC_INCLUDE} ${PETSC_CC_INCLUDES} -c -o  $@ $< ${SLEPC_EPS_LIB} 

${OBJ_DIR}/ex1.o: ${SRC_DIR}/ex1.c
	-${CC} ${SLEPC_INCLUDE} ${PETSC_CC_INCLUDES} -c -o $@ $< ${SLEPC_EPS_LIB}

${BIN_DIR}/ex1: ${OBJ_DIR}/read2.o ${LIB_DIR}/irpf90.a ${OBJ_DIR}/ex1.o ${SRC_DIR}/read2.h ${SRC_DIR}/stimsyr.h chkopts
	    -${CLINKER} ${SLEPC_INCLUDE} ${PETSC_CC_INCLUDES} -o ${BIN_DIR}/ex1 ${OBJ_DIR}/ex1.o ${OBJ_DIR}/read2.o ${LIB_DIR}/irpf90.a ${SLEPC_EPS_LIB}# -lifcore -lirc -lcomposerxe_gen_helpers_core_2.3
#    ${RM} ex1.o read2.o
