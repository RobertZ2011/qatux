CPP := clang++
CFLAGS := `pkg-config --cflags eigen3`
INCLUDE := /usr/include/eigen3/unsupported

all:

libqatux.a:
	$(CPP) src/general.cpp src/single_gate.cpp src/double_gate.cpp src/Qatux.cpp -o qatux.o $(CFLAGS) -I $(INCLUDE)
