#ifndef LOGGER_CUMC_H
#define LOGGER_CUMC_H

#include <stdexcept>
#include <string>
#include <iostream>
#include <fstream>
#include <functional>
#include <boost/ref.hpp>
#include <boost/bind.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/barrier.hpp>

#include <boost/log/common.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/attributes.hpp>
#include <boost/log/sinks.hpp>
#include <boost/log/sources/logger.hpp>
#include <boost/log/utility/record_ordering.hpp>

namespace logging = boost::log;
namespace attrs = boost::log::attributes;
namespace src = boost::log::sources;
namespace sinks = boost::log::sinks;
namespace expr = boost::log::expressions;
namespace keywords = boost::log::keywords;

using boost::shared_ptr;

enum
{
    LOG_RECORDS_TO_WRITE = 10,
    THREAD_COUNT = 2
};

BOOST_LOG_INLINE_GLOBAL_LOGGER_DEFAULT(test_lg, src::logger_mt)
//! This function is executed in multiple threads
void thread_fun(boost::barrier& bar)
{
    // Wait until all threads are created
    bar.wait();

    // Here we go. First, identify the thread.
    BOOST_LOG_SCOPED_THREAD_TAG("ThreadID", boost::this_thread::get_id());

    // Now, do some logging
    for (unsigned int i = 0; i < LOG_RECORDS_TO_WRITE; ++i)
    {
        BOOST_LOG(test_lg::get()) << "Log record " << i;
    } 
}
typedef sinks::text_ostream_backend backend_t;
typedef sinks::asynchronous_sink<backend_t, sinks::unbounded_ordering_queue<
        logging::attribute_value_ordering< unsigned int, std::less< unsigned int > > > > sink_t;

namespace Util
{
    class Logger
    {
        public:
            Logger(){}
            
            ~Logger()
            {
                // Wait until all action ends
                threads.join_all();
                
                // Flush all buffered records
                sink->stop();
                sink->flush();

            }
            void Initialize(const std::string& log_file_dir, const std::string& app_name)
            {   
                auto log_file_path = log_file_dir + "/" + app_name;
                std::time_t now_time = std::time({});
                tm* now_time_ptr = std::localtime(&now_time);
                char dt_string[100];

                std::strftime(dt_string, 90, "%Y%m%d_%H.%M.%S", now_time_ptr);
                std::stringstream ss;
                ss << log_file_path << "_" << dt_string << ".log";
                std::cout << ss.str() << std::endl;
                
                // Open a rotating text file
                shared_ptr< std::ostream > strm(new std::ofstream("test.log"));
                if (!strm->good())
                    throw std::runtime_error("Failed to open a text log file");

                // Create a text file sink
                typedef sinks::text_ostream_backend backend_t;
                typedef sinks::asynchronous_sink<
                    backend_t,
                    sinks::unbounded_ordering_queue<
                        logging::attribute_value_ordering< unsigned int, std::less< unsigned int > >
                    >
                > sink_t;
                sink = boost::make_shared<sink_t>(new sink_t(
                    boost::make_shared< backend_t >(),
                    // We'll apply record ordering to ensure that records from different threads go sequentially in the file
                    keywords::order = logging::make_attr_ordering("RecordID", std::less< unsigned int >()),
                    keywords::auto_flush = true                    
                    ));

                sink->locked_backend()->add_stream(strm);
                
                sink->set_formatter
                (
                    expr::format("%1%: [%2%] [%3%] - %4%")
                        % expr::attr< unsigned int >("RecordID")
                        % expr::attr< boost::posix_time::ptime >("TimeStamp")
                        % expr::attr< boost::thread::id >("ThreadID")
                        % expr::smessage                    
                );

                // Add it to the core
                logging::core::get()->add_sink(sink);

                // Add some attributes too
                logging::core::get()->add_global_attribute("TimeStamp", attrs::local_clock());
                logging::core::get()->add_global_attribute("RecordID", attrs::counter< unsigned int >());

                // Create logging threads
                boost::barrier bar(THREAD_COUNT);
                for (unsigned int i = 0; i < THREAD_COUNT; ++i)
                    threads.create_thread(boost::bind(&thread_fun, boost::ref(bar)));

       
                // Wait until all action ends
                //threads.join_all();

                // Flush all buffered records
                //sink->stop();
                //sink->flush();


                /* 

                
                logging::add_file_log
                (
                    keywords::file_name = ss.str(),
                    keywords::format = "[%TimeStamp%] [%ThreadID%] [%Severity%] %Message%",
                    keywords::auto_flush = true
                );

                logging::add_common_attributes();
                 */
                /* boost::shared_ptr<logging::core> core = logging::core::get();
                // Create a backend and initialize it with a stream
                boost::shared_ptr<logging::sinks::text_ostream_backend > backend = boost::make_shared<logging::sinks::text_ostream_backend >();
                backend->add_stream(boost::shared_ptr< std::ostream >(&std::clog, boost::null_deleter()));
                
                // Wrap it into the frontend and register in the core
                boost::shared_ptr<sink_t> sink(new sink_t(backend));
                core->add_sink(sink);

                // You can manage filtering and formatting through the sink interface
                //sink->set_filter(logging::expressions::attr<severity_level>("Severity") >= warning); */
                
            }
            void Debug(const std::string& msg)
            {
                BOOST_LOG(test_lg::get()) << msg;
                //BOOST_LOG_TRIVIAL(debug) << msg;
            }
            void Info(const std::string& msg)
            {
                BOOST_LOG(test_lg::get()) << msg;
                //BOOST_LOG_TRIVIAL(info) << msg;
            }
            void Error(const std::string& msg)
            {
                BOOST_LOG(test_lg::get()) << msg;
                //BOOST_LOG_TRIVIAL(error) << msg;
            }
        private:
            boost::shared_ptr<sink_t> sink;           
            boost::thread_group threads;
    };
};

#endif
