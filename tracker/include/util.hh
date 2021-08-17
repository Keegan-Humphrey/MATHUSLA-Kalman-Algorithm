#include <iostream>
#include <dirent.h>
#include "TString.h"

#ifndef UTIL_HH
#define UTIL_HH

namespace io {

bool is_directory(std::string name){
    return (opendir(name.c_str()) == NULL) ? false : true;
}

std::vector<TString> ProcessDirectory(std::string directory, std::string current_path)
{
   std::vector<TString> files = {};

   std::string path_to_dir = current_path + directory;
   std::string new_path = path_to_dir + '/';
 
   auto dir = opendir(path_to_dir.c_str());
   if (dir == NULL)
   {
      std::cout << "Could not open directory: " << directory.c_str() << std::endl;
      return {};
   }

   auto entity = readdir(dir);
   while (entity != NULL)
   {
      if(entity->d_type == DT_DIR) 
      {
          if(entity->d_name[0] == '.') {
            entity = readdir(dir);
            continue;
          }
            
          auto internal_file_names = ProcessDirectory(std::string(entity->d_name), new_path);
          for (auto file : internal_file_names) files.push_back(file);  

      } else if(entity->d_type == DT_REG){
          TString file_name = TString(entity->d_name);
          if (file_name.EndsWith(".root") and !(file_name.Contains("temp"))) files.push_back( TString(new_path) + TString(entity->d_name) );
          //ProcessFile(std::string(entity->d_name), new_path);
      }

    entity = readdir(dir);

      
   }
   
   closedir(dir);

   return files;

}


}; //namespace io

#endif