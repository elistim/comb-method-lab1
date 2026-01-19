CC = gcc
CFLAGS = -Iinclude -Wall
SRC = src/main.c
BUILD_DIR = build
OBJ = $(BUILD_DIR)/main.o
TARGET = $(BUILD_DIR)/my_program

SHELL = cmd.exe
export PATH := C:/msys64/mingw64/bin;$(PATH)

all: $(TARGET)

$(BUILD_DIR):
	cmd /c "if not exist $(BUILD_DIR) mkdir $(BUILD_DIR)"

$(TARGET): $(OBJ)
	$(CC) -o $@ $^

$(BUILD_DIR)/main.o: src/main.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	cmd /c "if exist $(BUILD_DIR) rd /S /Q $(BUILD_DIR)"

.PHONY: all clean