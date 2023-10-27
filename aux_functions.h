#ifndef AUX_FUNCTIONS_H
#define AUX_FUNCTIONS_H

#include <string>
#include <memory>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include <codecvt>
#include <chrono>

typedef std::unique_ptr<std::vector<std::string>> StrListPtr;
typedef std::unique_ptr<std::string> stringptr;

namespace Util
{
    const std::unordered_map<char, uint8_t> char_to_num_map {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};
    const std::unordered_map<char, uint8_t> complement_char_to_num_map {{'T', 0}, {'G', 1}, {'C', 2}, {'A', 3}};
    const char num_to_char_map[]{'A', 'C', 'G', 'T'};
    const std::unordered_map<char, char> char_complement_map {{'A', 'T'}, {'C', 'G'}, {'G', 'C'}, {'T', 'A'}};

    StrListPtr split(const std::string& str, char delim);

    int convert_to_int(const std::string& str);
    
    bool convert_to_bool(const std::string& str);

    std::string convert_to_string(int o);
    std::string convert_to_string(bool o);
    std::string convert_to_string(unsigned int o);
    std::string convert_to_string(u_int64_t o);
    void trim_str(std::string& str);
    void remove_first_chars(std::string& str, int num_of_chars);
    void remove_last_chars(std::string& str, int num_of_chars);
    void get_first_chars(std::string& str, int num_of_chars);
    void get_last_chars(std::string& str, int num_of_chars);
    uint8_t char_to_int(char, bool=false);
    //char num_to_char(uint8_t n);
    char complement(char c);
   // uint32_t pack(std::string::const_iterator begin, std::string::const_iterator end, bool complement);
    std::chrono::system_clock::time_point start_timing();
    std::chrono::microseconds stop_timing(std::chrono::system_clock::time_point start);

    template<typename T>
    std::string convert_to_csv(std::vector<T> v)
    {
        std::string result;
        for (T item: v)
        {
            result.append(convert_to_string(item));
            result.append(",");
        }

        auto last_comma = result.find_last_of(',');
        if (last_comma != std::string::npos)
            result.erase(last_comma);

        return result;
    }
}

#endif