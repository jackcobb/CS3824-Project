#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=Cygwin_4.x-Windows
CND_DLIB_EXT=dll
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/DnaSequenceVector.o \
	${OBJECTDIR}/IDnaRepository.o \
	${OBJECTDIR}/RandomEnumerator.o \
	${OBJECTDIR}/main.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/motiffinder.exe

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/motiffinder.exe: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/motiffinder ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/DnaSequenceVector.o: DnaSequenceVector.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I. -include IDnaSequence.h -include Neocleotide.h -include RandomEnumerator.h -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/DnaSequenceVector.o DnaSequenceVector.cpp

${OBJECTDIR}/IDnaRepository.o: IDnaRepository.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I. -include IDnaSequence.h -include Neocleotide.h -include RandomEnumerator.h -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/IDnaRepository.o IDnaRepository.cpp

${OBJECTDIR}/RandomEnumerator.o: RandomEnumerator.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I. -include IDnaSequence.h -include Neocleotide.h -include RandomEnumerator.h -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/RandomEnumerator.o RandomEnumerator.cpp

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I. -include IDnaSequence.h -include Neocleotide.h -include RandomEnumerator.h -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/motiffinder.exe

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
