#ifndef CONFIG_READER_H
#define CONFIG_READER_H
#include <string>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <memory>


namespace Util
{
    class ConfigReader
    {
        public:
            ConfigReader(){}
            void Initialize(const std::string& file_path);
            template<typename T>
            T GetValue(const std::string& key)
            {
                if (m_config_map.find(key) != m_config_map.end())
                {
                    std::stringstream ss{m_config_map[key]};
                    T t;
                    ss >> t;
                    return t;
                }
                else 
                    throw key + " does not exist!";
            }
            std::unique_ptr<std::string> GetAllKeyValues();

        private:
            
            std::unordered_map<std::string, std::string> m_config_map;
            
    };

    

}


#endif
