CC = C:/msys64/mingw64/bin/gcc
CFLAGS = -Iinclude -Wall
SRC = src/main.c
OBJ = $(SRC:.c=.o)
TARGET = my_program

SHELL = cmd.exe

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) -o $@ $^

src/main.o: src/main.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	cmd /c "if exist src\\main.o del /Q src\\main.o"
	cmd /c "if exist $(TARGET).exe del /Q $(TARGET).exe"