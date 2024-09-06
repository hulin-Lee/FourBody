#include "ArgumentParser.h"
#include <cstdlib>  // for std::atoi, std::atof
#include <iostream>

ProgramOptions parseArguments(int argc, char *argv[]) {
    ProgramOptions options;
    // 初始化默认值
    options.seed = -1;
    options.q_MBHB = -1.0;
    options.Mp = -1.0;
    options.Ms = -1.0;
    options.a_CM = -1.0;
    options.a_b = -1.0;
    options.I_CM = -1.0;
    options.tau = -1.0;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg.find("--seed=") == 0) {
            options.seed = std::atoi(arg.substr(7).c_str());
        } else if (arg.find("--M_IMBH=") == 0) {
            options.q_MBHB = std::atof(arg.substr(9).c_str());
        } else if (arg.find("--Mp=") == 0) {
            options.Mp = std::atof(arg.substr(5).c_str());
        } else if (arg.find("--Ms=") == 0) {
            options.Ms = std::atof(arg.substr(5).c_str());
        } else if (arg.find("--a_CM=") == 0) {
            options.a_CM = std::atof(arg.substr(7).c_str());
        } else if (arg.find("--a_b=") == 0) {
            options.a_b = std::atof(arg.substr(5).c_str());
        } else if (arg.find("--I_CM=") == 0) {
            options.I_CM = std::atof(arg.substr(7).c_str());
        } else if (arg.find("--tau=") == 0) {
            options.tau = std::atof(arg.substr(6).c_str());
        }
    }

    return options;
}
