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

int Util::convert_to_int(const std::string& str)
{
    std::stringstream ss(str);
    int num;
    ss >> num;
    return num;
}
    
bool Util::convert_to_bool(const std::string& str)
{
    std::stringstream ss(str);
    bool b;
    ss  >> std::boolalpha >> b;
    return b;
}

void Util::remove_first_chars(std::string& str, int num_of_chars)
{
    if (num_of_chars > 0)
        str.erase(0, num_of_chars);
}

void Util::remove_last_chars(std::string& str, int num_of_chars)
{
    auto remove_idx = str.size() - num_of_chars;
    if (remove_idx < str.size() && remove_idx >= 0)
        str.erase(remove_idx, num_of_chars);
}

void Util::get_first_chars(std::string& str, int num_of_chars)
{
    Util::remove_last_chars(str, str.size() - num_of_chars);
}

void Util::get_last_chars(std::string& str, int num_of_chars)
{
    Util::remove_first_chars(str, str.size() - num_of_chars);
}


std::string Util::convert_to_string(bool o)
{
    if (o) return "True";
    return "False";
}

 std::chrono::system_clock::time_point Util::start_timing()
 {
    return std::chrono::high_resolution_clock::now();
 }
 std::chrono::microseconds Util::stop_timing(std::chrono::system_clock::time_point start)
 {
    auto stop = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
 }