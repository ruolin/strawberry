/*
>HEADER
    Copyright (c) 2015 Ruolin Liu rliu0606@gmail.com
    This file is part of Strawberry.
    Strawberry is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Strawberry is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Strawberry.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
*/




#include "logger.hpp"
void file_log_policy::open_ostream(const std::string& name)
 {
  out_stream->open(name.c_str(), std::ios_base::binary|std::ios_base::out );
  if( !out_stream->is_open() )
  {
   throw(std::runtime_error("LOGGER: Unable to open an output stream"));
  }
 }

 void file_log_policy::close_ostream()
 {
  if( out_stream )
  {
   out_stream->close();
  }
 }

 void file_log_policy::write(const std::string& msg)
 {
  (*out_stream)<<msg<<std::endl;
 }

