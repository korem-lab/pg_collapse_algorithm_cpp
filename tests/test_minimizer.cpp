#include "../global_copan.h"
#include "../config_reader.h"
#include "../pe_ext_reader.h"
#include "../ext_contig.h"
#include "../async_logger.h"
#include "../aux_functions.h"
#include "../minimizer.h"
#include <string>
#include <iostream>
#include <fstream>
#include <ostream>

 Util::ConfigReader config;
 Util::AsyncLogger logger;

void TestKmerCount(ContigContainerPtr contigs)
{
       // count the kmers in the input
    logger.Info("counting kmers...");
    std::unordered_map<unsigned int, unsigned int> kc(contigs->size()* 1500);
    auto k = config.GetValue<int>("kmer_size");
    for (const ExtContig &c : *contigs)
    {
        mzr::count_kmers(c.Seq, kc, k);
    }
    auto test_size = kc.size();
    logger.Info("total kmers: " + Util::to_str(kc.size())); 

    std::string afile = "tests/data/kmers_py.txt";
    std::ifstream kmers_file(afile, std::ios::out);
    logger.Info("Reading kmers from file " + afile);
    long count = 0;
    long missing = 0;
    if(kmers_file.is_open())
    {
        while(kmers_file)
        {
            std::string line;
            std::getline(kmers_file, line);
            auto parts = Util::split(line, ' ');
            if (parts->size() == 2)
            {
                auto key = Util::to_int(parts->at(0));
                auto val = Util::to_int(parts->at(1));
                count++;
                if (kc.find(key) == kc.end())
                    missing++;
            }
            if (count % 10000000 == 0)
                logger.Info(" * Read " + Util::to_str(count) + " records.");
        }
    }
    kmers_file.close();   

    logger.Info("Read a total of " + Util::to_str(count) + " records.");

    if (test_size != count)
    {
        logger.Error("Total number of kmers do not match.");
        if (missing > 0)
        {
            logger.Error("Solution is missing " + Util::to_str(missing) + " records");
        }
    }
    else
    {
        logger.Info("Kmers match.");
    }
}
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
    
    TestKmerCount(ccptr);
    
    
    Util::Close();
    return 0;
}