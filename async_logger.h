#ifndef ASYNC_LOGGER_H
#define ASYNC_LOGGER_H


#include <iostream>
#include <thread>
#include <queue>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
#include <mutex>
#include "aux_functions.h"

namespace Util
{
namespace LogSeverity
{
    const std::string DebugLevel = "[debug]";
    const std::string InfoLevel = "[info]";
    const std::string ErrorLevel = "[error]";
};
class AsyncLogger {
private:
    std::queue<std::string> m_log_queue;
    std::thread m_logging_thread;
    std::ofstream m_log_file;
    bool m_stop_logging = false;
    std::mutex mutex;

    void enqueue(stringptr val) {
        std::lock_guard<std::mutex> lock(mutex);
        m_log_queue.push(*val);
    }

    stringptr dequeue() {
        std::lock_guard<std::mutex> lock(mutex);
        if (!m_log_queue.empty())
        {
            std::string val = m_log_queue.front();
            m_log_queue.pop();
            return std::make_unique<std::string>(val);
        }
        return std::make_unique<std::string>("");
    }

    std::string get_file_name(const std::string& log_file_dir, const std::string& app_name) 
    {
        auto log_file_path = log_file_dir + "/" + app_name;
        std::time_t now_time = std::time({});
        tm* now_time_ptr = std::localtime(&now_time);
        char dt_string[100];

        std::strftime(dt_string, 90, "%Y%m%d_%H.%M.%S", now_time_ptr);
        std::stringstream log_file_name;
        log_file_name << log_file_path << "_" << dt_string << ".log";
        return std::move(log_file_name.str());
    }
    
    stringptr format_message(const std::string& msg, const std::string& log_level)
    {
        std::time_t time = std::time(NULL);
        char buffer[256];
        std::strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S  ", std::localtime(&time));    
        std::stringstream ss;
        ss << buffer << std::this_thread::get_id() << "  " << log_level << "\t" << msg << std::endl; 
        return std::make_unique<std::string>(ss.str());
    }


public:
    void Debug(const std::string& message) {        
        enqueue(format_message(message, LogSeverity::DebugLevel));
    }

    void Info(const std::string& message){
        enqueue(format_message(message, LogSeverity::InfoLevel));
    }

    void Error(const std::string& message){
        enqueue(format_message(message, LogSeverity::ErrorLevel));
    }
    
    void Debug()
    {
        if (!m_log_queue.empty()) {                
                std::string message = *dequeue(); 
                if (message != "")
                {        
                    std::cout << message;
                }
            }
    }
    void Start(const std::string& log_file_dir, const std::string& app_name) 
    {
        m_log_file.open(get_file_name(log_file_dir, app_name));
        m_logging_thread = std::thread([this]() {
        while (!m_stop_logging) {
            //std::this_thread::sleep_for(std::chrono::milliseconds(100));
            if (!m_log_queue.empty()) {                
                std::string message = *dequeue(); 
                if (message != "")
                {        
                    std::cout << message;
                    m_log_file.write(message.c_str(), message.length());
                }
            }
        }
        });  
    }

  void Stop() {
    while(true)
    {
        std::this_thread::sleep_for(std::chrono::microseconds(5));
        if (m_log_queue.empty())
        {
            m_stop_logging = true;
            m_logging_thread.join();
            m_log_file.close();
            break;
        } 
    }
  }
};
};
#endif