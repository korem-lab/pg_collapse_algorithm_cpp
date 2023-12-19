#include "../global_copan.h"
#include "../config_reader.h"
#include "../pe_ext_reader.h"
#include "../ext_contig.h"
#include "../async_logger.h"
#include "../aux_functions.h"
#include <string>
#include <iostream>

 Util::ConfigReader config;
 Util::AsyncLogger logger;

int main(int argc, char *argv[])
{
    auto config_name = Util::ParseConfigName(argc, argv);
    //Read from config file and initialize a logger
    Util::InitializeServices(config_name);    
    
    PeExtReader per{};
    ContigContainerPtr ccptr(new ContigContainer());
    IdMapPtr imptr(new IdMap());
    // Reads fasta file and populates contig container and idmap
    auto num_samples = per.Parse(ccptr, imptr);

    logger.Info("Processed " + Util::to_str(num_samples) + " sample(s).");
    Util::Close();

    return 0;
}