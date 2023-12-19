#include "../global_copan.h"
#include "../config_reader.h"
#include "../pe_ext_reader.h"
#include "../ext_contig.h"
#include "../async_logger.h"
#include "../aux_functions.h"
#include "../alignment.h"
#include <string>
#include <iostream>
#include <ostream>

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
    auto alignments = compute_alignments(ccptr); 
    std::vector<PyAlignment>::const_iterator sol = alignments.cbegin();
    std::string afile = "tests/data/alignments_py.txt";
    std::ifstream pya_file(afile);
    std::vector<std::string> missing_alignments{};
    long count = 0;
    if (pya_file.is_open())
    {
        while(pya_file)
        {
            std::string line;
            std::getline(pya_file, line);
            auto parts = Util::split(line, ' ');

            if (parts->size() >= 6)
            {
                count++;
                PyAlignment p;
                p.q_begin = Util::to_int(parts->at(0));
                p.q_end = Util::to_int(parts->at(1));
                p.q_cid = Util::to_int(parts->at(2));
                p.s_begin = Util::to_int(parts->at(3));
                p.s_end = Util::to_int(parts->at(4));
                p.s_cid = Util::to_int(parts->at(5));

                if (sol != alignments.cend())
                {
                    if (!(sol->q_begin == p.q_begin && sol->q_end == p.q_end &&
                        sol->q_cid == p.q_cid && sol->s_begin == p.s_begin &&
                        sol->s_end == p.s_end && sol->s_cid == p.s_cid))
                        missing_alignments.push_back("MISMATCH: " + sol->ToString());                    
                    sol++;
                }
                else
                {
                    missing_alignments.push_back(p.ToString());
                }
            }

        }
    }
    pya_file.close();

    logger.Info("Finished reading alignments file " + afile);
    logger.Info("Read " + Util::to_str(count) + " records.");
    if (count != alignments.size())
    {        
        logger.Error("Alignments don't match");
        if (missing_alignments.size() > 0)
        {
            std::string bfile = "tests/output/missing_alignments.txt";
            std::ofstream al_file(bfile, std::ios::out);
            if (al_file.is_open())
            {
                for(const auto& a : missing_alignments)
                {
                    al_file << a << "\n";
                }
            }
            al_file.close();
        }
    }
    else
    {
        logger.Info("Alignments match.");
    }
    
    Util::Close();
    return 0;
}