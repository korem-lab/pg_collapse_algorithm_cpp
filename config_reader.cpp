#include "config_reader.h"
#include "aux_functions.h"


using namespace Util;

void split_config_line(const std::string& config_line, std::string& key, std::string& value, char delim) {

    auto parts = split(config_line, delim);
    if (parts->size() == 2)
    {
        key = parts->at(0);
        value = parts->at(1);
    }
    else if (parts->size() == 1)
    {
        key = parts->at(0);
    }	
}

void remove_comment(std::string& config_line)
{
    auto comment_pos = config_line.find('#');
    if (comment_pos != std::string::npos)
        config_line = config_line.substr(0, comment_pos);
}

void ConfigReader::Initialize(const std::string& file_path)
{
    try
    {
        std::ifstream configFile(file_path);
        if(configFile.is_open())
        {
            std::string line, key, value;
            char delim = '=';
            while(configFile)
            {
                std::getline(configFile, line);
                remove_comment(line);
                if (line.size() == 0) // skip comments
                    continue;
                
                split_config_line(line, key, value, delim);
                if (m_config_map.find(key) == m_config_map.end())
                {
                    m_config_map[key] = value;
                }
            }
            configFile.close();
        }
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';     

        throw e;
    }    
}

std::unique_ptr<std::string> ConfigReader::GetAllKeyValues()
{
    std::string s("\n>>>> config.ini <<<< \n");
    s.append("--------------------\n");
    for(auto keyVal: m_config_map)
    {
        s.append(keyVal.first + "=" + keyVal.second + "\n");
    }
    return std::make_unique<std::string>(s);
}

template<>
bool ConfigReader::GetValue<bool>(const std::string& key)
{
    if (m_config_map.find(key) != m_config_map.end())
    {
        std::stringstream ss{m_config_map[key]};
        bool t;
        ss >> std::boolalpha >> t;
        return t;
    }
    else 
        throw key + " does not exist!";
}