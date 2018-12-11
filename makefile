all: gtshark

GTShark_ROOT_DIR = .
GTShark_MAIN_DIR = src
LIBS_DIR = . #/usr/local/lib
INCLUDE_DIR= . #/usr/local/include
HTS_INCLUDE_DIR=htslib/include
HTS_LIB_DIR=htslib/lib

CC 	= g++
CFLAGS	= -Wall -O3 -m64 -std=c++11 -pthread -mavx -I $(HTS_INCLUDE_DIR) -I $(INCLUDE_DIR) -fpermissive
CLINK	= -lm -O3 -std=c++11 -pthread -mavx  -lhts -llzma -L $(HTS_LIB_DIR) -L $(LIBS_DIR)

ifdef MSVC     # Avoid the MingW/Cygwin sections
    uname_S := Windows
else                          # If uname not available => 'not' 
    uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
endif
ifeq ($(uname_S),Linux)
	CLINK+=-fabi-version=6
endif

# default install location (binary placed in the /bin folder)
prefix      = /usr/local

# optional install location
exec_prefix = $(prefix)


%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

gtshark: $(GTShark_MAIN_DIR)/application.o \
	$(GTShark_MAIN_DIR)/cfile.o \
	$(GTShark_MAIN_DIR)/lzma_wrapper.o \
	$(GTShark_MAIN_DIR)/main.o \
	$(GTShark_MAIN_DIR)/pbwt.o \
	$(GTShark_MAIN_DIR)/sfile.o \
	$(GTShark_MAIN_DIR)/utils.o \
	$(GTShark_MAIN_DIR)/vcf.o 
	$(CC) -o $(GTShark_ROOT_DIR)/$@  \
	$(GTShark_MAIN_DIR)/application.o \
	$(GTShark_MAIN_DIR)/cfile.o \
	$(GTShark_MAIN_DIR)/lzma_wrapper.o \
	$(GTShark_MAIN_DIR)/main.o \
	$(GTShark_MAIN_DIR)/pbwt.o \
	$(GTShark_MAIN_DIR)/sfile.o \
	$(GTShark_MAIN_DIR)/utils.o \
	$(GTShark_MAIN_DIR)/vcf.o \
	$(CLINK)

clean:
	-rm $(GTShark_MAIN_DIR)/*.o
	-rm gtshark

install:
	mkdir -p -m 755 $(exec_prefix)/bin
	cp gtshark $(exec_prefix)/bin/
	
uninstall:
	rm  $(exec_prefix)/bin/gtshark
	

