// UtilStats.h
#pragma once

struct UtilStats {
    double total_core = 0.0;      // Σ sum_core
    double total_capacity = 0.0;  // Σ (P * T_wall)
    long long calls = 0;

    inline void add(double sum_core, double T_wall, int P) {
        total_core     += sum_core;
        total_capacity += T_wall * (double)P;
        calls++;
    }

    inline double utilization() const {
        return total_capacity > 0.0 ? total_core / total_capacity : 0.0;
    }
};
