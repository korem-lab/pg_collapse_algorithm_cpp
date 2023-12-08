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
    //const char num_to_char_map[]{'A', 'C', 'G', 'T'};
    
    StrListPtr split(const std::string& str, char delim);

    int convert_to_int(const std::string& str);
    
    bool convert_to_bool(const std::string& str);
     
    template<typename T>
    std::string convert_to_string(T o)
    {
        std::stringstream ss;
        ss << o;
        return ss.str();
    }  
    std::string convert_to_string(bool o);

    void trim_str(std::string& str);
    void remove_first_chars(std::string& str, int num_of_chars);
    void remove_last_chars(std::string& str, int num_of_chars);
    void get_first_chars(std::string& str, int num_of_chars);
    void get_last_chars(std::string& str, int num_of_chars);  

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
    template<typename T>   
    std::string convert_to_csv(std::unordered_map<T, char> map)
    {
        std::string result;
        for(const auto& [key, val]: map)
        {
            result.append("{");
            result.append(convert_to_string(key));
            result.append(": '");
            result.append(convert_to_string(val));
            result.append("'},");
        }

        auto last_comma = result.find_last_of(',');
        if (last_comma != std::string::npos)
            result.erase(last_comma);

        return result;
    }     
    
}

#endif