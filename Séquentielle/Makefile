CXX := mpicxx
CXXFLAGS := -O3 -g
LDFLAGS := 

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
