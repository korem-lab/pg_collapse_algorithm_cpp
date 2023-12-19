#include "aux_functions.h"

StrListPtr Util::split(const std::string& str, char delim)
{
    std::stringstream ss(str);
    std::string word;
    std::vector<std::string> word_list;

    while (!ss.eof()) 
    {
        std::getline(ss, word, delim);
        Util::trim_str(word);
        word_list.push_back(word);
    }   
    return std::make_unique<std::vector<std::string>>(word_list);
}

int Util::to_int(const std::string& str)
{
    std::stringstream ss(str);
    int num;
    ss >> num;
    return num;
}
    
bool Util::to_bool(const std::string& str)
{
    std::stringstream ss(str);
    bool b;
    ss  >> std::boolalpha >> b;
    return b;
}