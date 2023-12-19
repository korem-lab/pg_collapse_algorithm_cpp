#include "global_copan.h"



void Util::InitializeServices(const char* config_name)
{
    try
    {
        config.Initialize(config_name);

        auto log_file_dir = config.GetValue<std::string>("log_file_dir");
        auto app_name = config.GetValue<std::string>("app_name");

        logger.Start(log_file_dir, app_name);
        logger.Info("Initializing... ");
        logger.Info(*config.GetAllKeyValues());
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        logger.Stop();
    }
}

const char* Util::ParseConfigName(int argc, char *argv[])
{
    if (argc > 1)
        return argv[1];
    return "config.ini";
}

void Util::Close()
{
    logger.Info("Complete!"); 
    logger.Stop();
}