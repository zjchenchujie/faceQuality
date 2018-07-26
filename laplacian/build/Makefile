##
## ****************************************************************************
##
##            Copyright(c), HIKDATA.All Rights Reserved
##
## file         : 
##              :
## brief        : 
##              :
##              :
## author       :duxin@hikdata.com
## version      :
## change history:
## ****************************************************************************
##


## ----------------------------------------------------------------------------
## module name
## ----------------------------------------------------------------------------

## 指定平台
PLATFORM=linux64_gcc

MODULE_NAME = test

## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------
ROOT_DIR = ..

## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------
include $(ROOT_DIR)/build/platform.rules


##
SRC_DIR = $(ROOT_DIR)/src

##

INC_DIR = $(ROOT_DIR)/include
#INC_DIR += -I/usr/include/python2.7
#INC_DIR += -I/usr/lib/python2.7/dist-packages/numpy/core/include
#INC_DIR += -I/usr/local/include
#INC_DIR += -I/usr/include/hdf5/serial
#INC_DIR += -I/usr/local/cuda/include

#INC_DIR += $(ROOT_DIR)/.build_release/src
#INC_DIR += /usr/local/cuda

##
TMP_DIR = ./tmp

##
# LIB_DIR = $(ROOT_DIR)/lib
LIB_DIR_OPENCV = /usr/local/lib/

## 
LDFLAGS  = -L$(LIB_DIR)
#LDFLAGS += -L$(LIB_DIR) -L/usr/lib -L/usr/local/lib -L/usr/lib -L/usr/lib/x86_64-linux-gnu/hdf5/serial -L/usr/local/cuda/lib64
#LDFLAGS += -L/usr/local/cuda/lib

#LDFLAGS += -L$(LIB_DIR) -l_face_align
#LDFLAGS += -L$(LIB_DIR) -l_aia_base

##
#LDFLAGS += -L$(LIB_DIR) -Wl,--whole-archive -l_face_recognise -Wl,--no-whole-archive

##
#LDFLAGS += -lcudart -lcublas -lcurand -lglog -lgflags -lprotobuf -lboost_system -lboost_filesystem -lm -lhdf5_hl -lhdf5 -lleveldb -lsnappy -llmdb
#LDFLAGS += -lopencv_imgcodecs -lboost_thread -lstdc++ -lcudnn -lboost_python -lpython2.7 -lopenblas
LDFLAGS += -lopencv_highgui
LDFLAGS += -lopencv_imgproc
LDFLAGS += -lopencv_core
LDFLAGS += -lopencv_imgcodecs
#LDFLAGS += -lcaffe
#LDFLAGS += -l_track_mot
##
IFLAGS += -I$(INC_DIR)

#CFLAGS +=-DHAVE_OPENCV


## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------
SRCS = $(SRC_DIR)/laplacian_main.cpp
SRCS += $(SRC_DIR)/laplacian.cpp
# SRCS += $(SRC_DIR)/hello.cpp
#SRCS += $(SRC_DIR)/aia_det.cpp
#SRCS += $(SRC_DIR)/aia_yolo_v2_output.cpp
       
OBJS = $(patsubst %.cpp,%.o,$(patsubst %.c,%.o,$(SRCS)))
OBJT = $(addprefix $(TMP_DIR)/,$(notdir $(OBJS)))


## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------
TARGET = $(ROOT_DIR)/bin/$(MODULE_NAME)

##	
all: tmpdir $(TARGET)

##
tmpdir:
	mkdir -p $(TMP_DIR)
	echo $(SRCS)


## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------
#.c.o:
#	$(CC) $(CFLAGS) $(IFLAGS) -c $^ -o $(addprefix $(TMP_DIR)/,$(@F)) $(OUTPUT)

## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------
.cpp.o:
	$(CPP) $(CFLAGS) $(IFLAGS) -c $^ -o $(addprefix $(TMP_DIR)/,$(@F)) $(OUTPUT)
	
$(TARGET):$(OBJS)
	$(CPP) -o $@ $(OBJT) $(LDFLAGS) -I$(INC_DIR)
	


## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------
clean:
	find . -type f \( -name 'core'                                              \
	               -o -name '*.bak'                                             \
	               -o -name '*~'                                                \
	               -o -name '*.o'                                               \
	               -o -name '.depend'                                           \
	               -o -name 'build.output'                                      \
	               \) -print | xargs rm -f
	rm -f $(OBJS)
	

	


## ----------------------------------------------------------------------------
## ---------------------------------------------------------------------------
.PHONY: all clean $(TARGET)



