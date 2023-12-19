#ifndef AUX_FUNCTIONS_H
#define AUX_FUNCTIONS_H

#include <string>
#include <memory>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <codecvt>
#include <chrono>
#include <sys/stat.h>
#include <algorithm> 

typedef std::unique_ptr<std::vector<std::string>> StrListPtr;
typedef std::unique_ptr<std::string> stringptr;

namespace Util
{  
    StrListPtr split(const std::string& str, char delim);

    int to_int(const std::string& str);
    
    bool to_bool(const std::string& str);     
    
    inline std::string to_str(bool o)
    {
        if (o) return "True";
        return "False";
    }
    
    inline void remove_first_chars(std::string& str, int num_of_chars)
    {
        if (num_of_chars > 0)
            str.erase(0, num_of_chars);
    }
    inline void remove_last_chars(std::string& str, int num_of_chars)
    {
        auto remove_idx = str.size() - num_of_chars;
        if (remove_idx < str.size() && remove_idx >= 0)
            str.erase(remove_idx, num_of_chars);
    }
    inline void get_first_chars(std::string& str, int num_of_chars)
    {
        Util::remove_last_chars(str, str.size() - num_of_chars);
    }

    inline void get_last_chars(std::string& str, int num_of_chars)
    {
        Util::remove_first_chars(str, str.size() - num_of_chars);
    }
    inline std::chrono::system_clock::time_point start_timing()
    {
        return std::chrono::high_resolution_clock::now();
    }
    inline std::chrono::microseconds stop_timing(std::chrono::system_clock::time_point start)
    {
        auto stop = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    }

    inline bool file_exists(const std::string& path)
    {
        struct stat buffer;
        return (stat(path.c_str(), &buffer) == 0);
    }
    inline std::string get_file_ext(const std::string& filename)
    {
        return filename.substr(filename.find_last_of('.') + 1);
    }
    
    // trim from start (in place)
    inline void ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
            return !std::isspace(ch);
        }));
    }

    // trim from end (in place)
    inline void rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
            return !std::isspace(ch);
        }).base(), s.end());
    }

    inline void trim_str(std::string& str)
    {
        rtrim(str);
        ltrim(str);
    }

        template<typename T>
    std::string convert_to_csv(std::vector<T> v)
    {
        std::string result;
        for (T item: v)
        {
            result.append(to_str(item));
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
            result.append(to_str(key));
            result.append(": '");
            result.append(to_str(val));
            result.append("'},");
        }

        auto last_comma = result.find_last_of(',');
        if (last_comma != std::string::npos)
            result.erase(last_comma);

        return result;
    } 
    template<typename T>
    std::string to_str(T o)
    {
        std::stringstream ss;
        ss << o;
        return ss.str();
    }  
}

#endif