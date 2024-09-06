# 定义编译器
CXX = g++

# 编译标志
CXXFLAGS = -std=c++11 -Wall -Wno-reorder -Wno-unused-but-set-variable -Wno-sign-compare -Wno-unused-variable

# 目标文件
TARGET = FourBody

# 源文件
SRCS = FourBody.cpp ArgumentParser.cpp

# 头文件
HEADERS = ArgumentParser.h nr3.h odeint.h funs.h stepperdopr853.h integration.h

# 默认规则：生成可执行文件
all: $(TARGET)

# 链接和生成可执行文件的规则
$(TARGET): $(SRCS) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRCS)

# 清理生成的文件
clean:
	rm -f $(TARGET)
