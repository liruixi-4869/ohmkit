CC      = gcc
CFLAGS  = -Wall -Wextra -O2
LDFLAGS = -lm
TARGET  = ohmkit
SRC     = ohmkit.c

.PHONY: all clean run test

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(TARGET)

run: $(TARGET)
	./$(TARGET)

test: $(TARGET)
	@echo "=== 表达式 ==="
	./$(TARGET) "4+((9||27)+(1||27))||27"
	@echo ""
	@echo "=== 并联 ==="
	./$(TARGET) "10 || 20 || 30"
	@echo ""
	@echo "=== Δ→Y ==="
	./$(TARGET) delta 10 20 30
	@echo ""
	@echo "=== 桥式电路 ==="
	./$(TARGET) bridge 1 2 3 4 5
	@echo ""
	@echo "=== 色环解码 ==="
	./$(TARGET) color red violet yellow gold
	@echo ""
	@echo "=== 色环编码 ==="
	./$(TARGET) findcolor 4700
