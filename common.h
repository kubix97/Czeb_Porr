#pragma once

constexpr int MEMALIGN = 16;

// 0 - Sequential mode, 1 - use SSE, 2 - use OpenMP
#define EXMODE	1 