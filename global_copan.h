#ifndef GLOBAL_COPAN_H
#define GLOBAL_COPAN_H

#include "config_reader.h"
#include "async_logger.h"

extern Util::ConfigReader config;
extern Util::AsyncLogger logger;

namespace Util
{
    void InitializeServices(const char* config_name);
    const char* ParseConfigName(int argc, char *argv[]);
    void Close();
};

#endif