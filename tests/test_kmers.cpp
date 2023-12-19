#include "../minimizer.h"
#include "../aux_functions.h"
#include <iostream>
#include <iomanip>

 Util::ConfigReader config;
 Util::AsyncLogger logger;
 
void test_aux_functions_get_first_chars()
{
    std::string str1{"AATTCCGGTT"};
    Util::get_first_chars(str1, 2);
    if (str1 == "AA")
        std::cout << "test1 of aux_functions::get_first_chars passed" << std::endl;
    else
        std::cout << "test1 of aux_functions::get_first_chars failed: " + str1 << std::endl;
    
    std::string str2{"AATTCCGGTT"};
    Util::get_first_chars(str2, str2.size());
    if (str2 == "AATTCCGGTT")
        std::cout << "test2 of aux_functions::get_first_chars passed" << std::endl;
    else
        std::cout << "test2 of aux_functions::get_first_chars failed" + str2 << std::endl;
}

void test_aux_functions_get_last_chars()
{
    std::string str1{"AATTCCGGTT"};
    Util::get_last_chars(str1, 2);
    if (str1 == "TT")
        std::cout << "test1 of aux_functions::get_last_chars passed" << std::endl;
    else
        std::cout << "test1 of aux_functions::get_last_chars failed: " + str1 << std::endl;
    
    std::string str2{"AATTCCGGTT"};
    Util::get_last_chars(str2, str2.size());
    if (str2 == "AATTCCGGTT")
        std::cout << "test2 of aux_functions::get_last_chars passed" << std::endl;
    else
        std::cout << "test2 of aux_functions::get_last_chars failed" + str2 << std::endl;
}

void test_aux_functions_remove_first_chars()
{
    std::string str1{"AATTCCGGTT"};
    Util::remove_first_chars(str1, 2);
    if (str1 == "TTCCGGTT")
        std::cout << "test1 of aux_functions::remove_first_chars passed" << std::endl;
    else
        std::cout << "test1 of aux_functions::remove_first_chars failed: " + str1 << std::endl;
    
    std::string str2{"AATTCCGGTT"};
    Util::remove_first_chars(str2, str2.size());
    if (str2 == "")
        std::cout << "test2 of aux_functions::remove_first_chars passed" << std::endl;
    else
        std::cout << "test2 of aux_functions::remove_first_chars failed: " + str2 << std::endl;
}

void test_aux_functions_remove_last_chars()
{
    std::string str1{"AATTCCGGTT"};
    Util::remove_last_chars(str1, 2);
    if (str1 == "AATTCCGG")
        std::cout << "test1 of aux_functions::remove_last_chars passed" << std::endl;
    else
        std::cout << "test1 of aux_functions::remove_last_chars failed: " + str1 << std::endl;
    
    std::string str2{"AATTCCGGTT"};
    Util::remove_last_chars(str2, str2.size());
    if (str2 == "")
        std::cout << "test2 of aux_functions::remove_last_chars passed" << std::endl;
    else
        std::cout << "test2 of aux_functions::remove_last_chars failed: " + str2 << std::endl;
} 

void test_kmer_pack()
{
    auto start = Util::start_timing();
    std::string s = "AATTCCGGTTAGCTGAATTCCGGTTAGCTGAATTCCGGTTAGCTGAATTCCGGTTAGCTG";
    int k = 15;
    mzr::KmerGenerator kmer_gen(s, k, 0x3FFFFFFF, true);
    auto km = kmer_gen.get_kmer();
    
    bool str_pack_is_correct = km.seq == 64404638;
    //bool rev_str_is_correct = kmer_gen.rev_seq == "GTCGATTGGCCTTAAGTCGATTGGCCTTAAGTCGATTGGCCTTAAGTCGATTGGCCTTAA";
    bool rev_str_pack_is_correct = kmer_gen.rev_kmer.seq == 309352975;
    

    while (!kmer_gen.empty()) 
        kmer_gen.get_kmer();

    auto min_kmer = kmer_gen.min_kmer_in_window(15);
    auto bool_is_min_kmer_seq_correct = min_kmer.seq == 64404638;

    auto dur = Util::stop_timing(start);

    if (str_pack_is_correct == rev_str_pack_is_correct == bool_is_min_kmer_seq_correct == true)
        std::cout << "Passed the test_kmer_pack test: " <<  dur.count() << std::endl;
    else    
        std::cout << "Failed the test_kmer_pack test: " <<  dur.count() << std::endl;
}

void test_min_kmer_in_window()
{
    std::string s = "AATTCCGGTTAGCTGAATTCCGGTTAGCTGAATTCCGGTTAGCTGAATTCCGGTTAGCTG";
    int k = 3;
    mzr::KmerGenerator kmer_gen(s, k, 0x3FFFFFFF, true);
    auto kmer1 = kmer_gen.get_kmer();
    uint8_t window = 1;
    auto km = kmer_gen.min_kmer_in_window(window);
    bool result = (km.seq == 3);
    if (result)
        std:: cout << "Passed the test_min_kmer_in_window " << std::endl;
    else
        std:: cout << "Failed the test_min_kmer_in_window " << std::endl;
}

int main()
{
    test_aux_functions_get_first_chars();
    test_aux_functions_get_last_chars();
    test_aux_functions_remove_first_chars();
    test_aux_functions_remove_last_chars();
    test_kmer_pack();
    test_min_kmer_in_window();
}

