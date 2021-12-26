#ifndef TEST_TOOLS_H
#define TEST_TOOLS_H


#include <string>
#include "gtest/gtest.h"

namespace arguments {
std::string testInputFile;
}

class testEnvironment : public testing::Environment {
 public:
  explicit testEnvironment(const std::string &command_line_arg) {
    arguments::testInputFile = command_line_arg;
  }
};

#endif