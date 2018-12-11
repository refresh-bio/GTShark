#pragma once
// *******************************************************************************************
// This file is a part of GTShark software distributed under GNU GPL 3 licence.
// The homepage of the GTShark project is https://github.com/refresh-bio/GTShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.0
// Date   : 2018-12-10
// *******************************************************************************************

#include <array>
#include <cstdint>

#define OUR_STRTOL

using namespace std;

typedef array<uint32_t, 4> stats_t;
typedef array<uint64_t, 4> stats64_t;

typedef uint64_t context_t;

typedef array<pair<uint8_t, uint32_t>, 2> run_t;

const uint32_t SIGMA = 4u;

// EOF
