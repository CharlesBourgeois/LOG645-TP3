CXX := mpicxx
CXXFLAGS := -O3 -g `pkg-config --cflags sfml-graphics sfml-window sfml-system`
LDFLAGS := `pkg-config --libs sfml-graphics sfml-window sfml-system`

TARGET := MPI

SRCS := MPI.cpp

OBJS := $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) -o $@ -c $< $(CXXFLAGS)

clean:
	rm -f $(TARGET) $(OBJS)

.PHONY: clean
