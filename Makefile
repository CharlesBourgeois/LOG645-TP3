CXX := mpicxx
CXXFLAGS := -O3 -g -I/usr/include  
LDFLAGS := -L/usr/lib -lsfml-graphics -lsfml-window -lsfml-system  

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
