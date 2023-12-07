CC := mpicc
CFLAGS := -O3 -g

TARGET := MPI
SRCS := main.c  

OBJS := $(SRCS:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

%.o: %.c
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -f $(TARGET) $(OBJS)

.PHONY: clean
