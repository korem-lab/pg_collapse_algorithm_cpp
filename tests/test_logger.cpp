#include "../async_logger.h"
#include "../global_copan.h"

 Util::ConfigReader config;
 Util::AsyncLogger logger;

int main(int argc, char *argv[])
{
    logger.Start("logs", "pg_collapse_algo");
    logger.Debug("Hello, world!");
    logger.Info("Hello, world!");
    logger.Error("Hello, world!");
    logger.Stop();

    return 0;
}  