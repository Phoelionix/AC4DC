/**
 * @file config_post_script.h
 * @brief 
 */
#include <string>
#pragma once

const std::string post_script = "~/AC4DC/config/post_scripts.sh"; // script to call after sim finished. Relative to binary (compiled program) path. (If "" does nothing.)
const bool execute_post_script = true;