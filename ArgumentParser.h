#ifndef ARGUMENTPARSER_H
#define ARGUMENTPARSER_H

#include <string>

// 参数结构体，包含需要解析的参数
struct ProgramOptions {
    int seed;
    double q_MBHB;
    double Mp;
    double Ms;
    double a_CM;
    double a_b;
    double I_CM;
    double tau;
};

// 参数解析函数的声明
ProgramOptions parseArguments(int argc, char *argv[]);

#endif // ARGUMENTPARSER_H
