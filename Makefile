include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
#CC=gcc
#FC = ifort
CLINKER = mpicc -fPIC  #-wd1572 -O3 -axAVX,SSE4.2 -fno-alias -no-prec-div -no-prec-sqrt -ip
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

${OBJ_DIR}/get_ntot.o: ${SRC_DIR}/get_ntot.c directories chkopts
	${CC} ${SLEPC_INCLUDE} ${PETSC_CC_INCLUDES} -c -o  $@ $< ${SLEPC_EPS_LIB} 

${OBJ_DIR}/read2.o: ${SRC_DIR}/read2.c directories chkopts
	${CC} ${SLEPC_INCLUDE} ${PETSC_CC_INCLUDES} -c -o  $@ $< ${SLEPC_EPS_LIB} 

${OBJ_DIR}/get_s2.o: ${SRC_DIR}/get_s2.c directories chkopts
	${CC} ${SLEPC_INCLUDE} ${PETSC_CC_INCLUDES} -c -o  $@ $< ${SLEPC_EPS_LIB} 

${OBJ_DIR}/get_s2_cyclic.o: ${SRC_DIR}/get_s2_cyclic.c directories chkopts
	${CC} ${SLEPC_INCLUDE} ${PETSC_CC_INCLUDES} -c -o  $@ $< ${SLEPC_EPS_LIB} 

${OBJ_DIR}/get_s2_mov.o: ${SRC_DIR}/get_s2_mov.c directories chkopts
	${CC} ${SLEPC_INCLUDE} ${PETSC_CC_INCLUDES} -c -o  $@ $< ${SLEPC_EPS_LIB} 

${OBJ_DIR}/get_dmat.o: ${SRC_DIR}/get_dmat.c directories chkopts
	${CC} ${SLEPC_INCLUDE} ${PETSC_CC_INCLUDES} -c -o  $@ $< ${SLEPC_EPS_LIB} 

${OBJ_DIR}/get_val_iaa2.o: ${SRC_DIR}/get_val_iaa2.c directories chkopts
	${CC} ${SLEPC_INCLUDE} ${PETSC_CC_INCLUDES} -c -o  $@ $< ${SLEPC_EPS_LIB} 

${OBJ_DIR}/ex1.o: ${SRC_DIR}/ex1.c
	-${CC} ${SLEPC_INCLUDE} ${PETSC_CC_INCLUDES} -c -o $@ $< ${SLEPC_EPS_LIB}

${BIN_DIR}/ex1: ${OBJ_DIR}/get_ntot.o ${OBJ_DIR}/read2.o ${OBJ_DIR}/get_s2_mov.o ${OBJ_DIR}/get_s2_cyclic.o ${OBJ_DIR}/get_s2.o ${OBJ_DIR}/get_dmat.o ${OBJ_DIR}/get_val_iaa2.o ${LIB_DIR}/irpf90.a ${OBJ_DIR}/ex1.o ${SRC_DIR}/read2.h ${SRC_DIR}/get_ntot.h ${SRC_DIR}/stimsyr.h chkopts
	    -${CLINKER} ${SLEPC_INCLUDE} ${PETSC_CC_INCLUDES} -o ${BIN_DIR}/ex1 ${OBJ_DIR}/ex1.o ${OBJ_DIR}/read2.o ${OBJ_DIR}/get_ntot.o ${OBJ_DIR}/get_s2.o ${OBJ_DIR}/get_s2_mov.o ${OBJ_DIR}/get_s2_cyclic.o ${OBJ_DIR}/get_dmat.o ${OBJ_DIR}/get_val_iaa2.o ${LIB_DIR}/irpf90.a ${SLEPC_EPS_LIB}
#    ${RM} ex1.o read2.o
