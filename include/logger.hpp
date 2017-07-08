
#ifndef STRAWBERRY_LOGGER_HPP_
#define STRAWBERRY_LOGGER_HPP_

#include "easylogging++.h"


inline void CreateLogger(const std::string& fname) {
   el::Configurations defaultConf;
   defaultConf.setToDefault();
   // Values are always std::string
   defaultConf.set(el::Level::Info,
           el::ConfigurationType::Format, "%datetime %level %msg");
   // default logger uses default configurations
   el::Loggers::reconfigureLogger("default", defaultConf);
   std::cerr << "Log to file " <<fname<<std::endl;
   // To set GLOBAL configurations you may use
   defaultConf.setGlobally(el::ConfigurationType::Filename, fname);
   defaultConf.setGlobally(el::ConfigurationType::ToStandardOutput, "false");
   defaultConf.setGlobally(el::ConfigurationType::ToFile, "true");
   el::Loggers::reconfigureLogger("default", defaultConf);
}

#endif /* STRAWBERRY_LOGGER_HPP_ */
