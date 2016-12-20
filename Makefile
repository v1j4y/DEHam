include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
CC=gcc-5
MAKE = /usr/bin/make
MKDIR_P = mkdir -p
OBJ_DIR = obj

irpf90.a:
	cd src && irpf90 init && $(MAKE) && cp irpf90.a ../

${OBJ_DIR}:
	${MKDIR_P} ${OBJ_DIR}

directories: ${OBJ_DIR}

obj/read2.o: src/read2.c directories
	${CC} -c -o $@ $<

obj/ex1.o: src/ex1.c
	-${CC} -c -o $@ $< ${SLEPC_EPS_LIB}

ex1: obj/read2.o irpf90.a obj/ex1.o src/read2.h src/stimsyr.h chkopts
	    -${CLINKER} -o ex1 obj/ex1.o obj/read2.o irpf90.a ${SLEPC_EPS_LIB}# -lifcore -lirc -lcomposerxe_gen_helpers_core_2.3
#    ${RM} ex1.o read2.o
