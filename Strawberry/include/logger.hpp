
#ifndef LOGGER_HPP_
#define LOGGER_HPP_

#include<memory>
#include<sstream>
#include<fstream>
#include <mutex>
//using namespace std;

/*
 * This file is modefied based on Dr.Dobb's lightweight c++ logger
 * http://www.drdobbs.com/cpp/a-lightweight-logger-for-c/240147505
 */

enum class severity_type
{
 debug = 1,
 error,
 warning,
 input_error
};

 class log_policy_interface
 {
 public:
  virtual void  open_ostream(const std::string& name) = 0;
  virtual void  close_ostream() = 0;
  virtual void  write(const std::string& msg) = 0;
  virtual ~log_policy_interface(){};
 };

 /*
  * Implementation which allow to write into a file
  */

 class file_log_policy : public log_policy_interface
 {
  std::unique_ptr< std::ofstream > out_stream;
 public:
  file_log_policy() : out_stream( new std::ofstream ) {}
  void open_ostream(const std::string& name);
  void close_ostream();
  void write(const std::string& msg);
  ~file_log_policy(){
     if( out_stream )
     {
      close_ostream();
     }
    }
 };

 template< typename log_policy >
 class logger
 {
  unsigned log_line_number;
  std::string get_time();
  std::string get_logline_header();
  std::stringstream log_stream;
  log_policy* policy;
  std::mutex write_mutex;

  //Core printing functionality
  void print_impl();
  template<typename First, typename...Rest>
  void print_impl(First parm1, Rest...parm);
 public:

  logger( const std::string& name );
  template< severity_type severity , typename...Args >
  void print( Args...args );
  ~logger();
 };

  template< typename log_policy >
  inline logger<log_policy>::logger( const std::string& name ) {
      policy = new log_policy();
      policy->open_ostream(name);
   }

 template< typename log_policy >
 inline logger< log_policy >::~logger()
 {
      policy->close_ostream();
      delete policy;
 }




 template< typename log_policy >
 inline void logger< log_policy >::print_impl()
 {
  policy->write( get_logline_header() + log_stream.str() );
  log_stream.str("");
 }

 template< typename log_policy >
  template<typename First, typename...Rest >
inline void logger< log_policy >::print_impl(First parm1, Rest...parm)
 {
  log_stream<<parm1;
  print_impl(parm...);
 }

 template< typename log_policy >
  template< severity_type severity , typename...Args >
inline  void logger< log_policy >::print( Args...args )
 {
  write_mutex.lock();
  switch( severity )
  {
   case severity_type::debug:
      log_stream<<"<DEBUG>: ";
      break;
   case severity_type::warning:
      log_stream<<"<WARNING>: ";
      break;
   case severity_type::error:
      log_stream<<"<ERROR>: ";
      break;
   case severity_type::input_error:
      log_stream<<"<INPUT_ERROR>: ";
      break;
  };
  print_impl( args... );
  write_mutex.unlock();
 }


 template< typename log_policy >
inline std::string logger< log_policy >::get_time()
 {
  std::string time_str;
  time_t raw_time;

  time( & raw_time );
  time_str = ctime( &raw_time );

  //without the newline character
  return time_str.substr( 0 , time_str.size() - 1 );
 }

 template< typename log_policy >
inline std::string logger< log_policy >::get_logline_header()
 {
  std::stringstream header;

  header.str("");
  header.fill('0');
  header.width(7);
  header << log_line_number++ <<" < "<<get_time()<<"> ~ ";
  return header.str();
 }


 static logger< file_log_policy > log_inst( "execution.log" );

 #define LOG log_inst.print< severity_type::debug >
 #define LOG_ERR log_inst.print< severity_type::error >
 #define LOG_WARN log_inst.print< severity_type::warning >
 #define LOG_INPUT log_inst.print<severity_type::input_error>

//static logger<file_log_policy> log_track("tracking.log");
//#define LOG_TRACK log_track.print< severity_type::debug >

#endif /* LOGGER_HPP_ */
