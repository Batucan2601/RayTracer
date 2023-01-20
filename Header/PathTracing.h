#pragma once
#include "../System_Files/support_files/parser.h"
#include <random>

// generator for importance sampling

// epsilon1 
std::random_device rd_importance_sampling_eps1;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen_importance_eps1(rd_importance_sampling_eps1()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis_importance_sampling_eps1(0 ,1);

//epslion2 
std::random_device rd_importance_sampling_eps2;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen_importance_eps2(rd_importance_sampling_eps2()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis_importance_sampling_eps2(0 ,1);

//generator for russian roulette 
std::random_device rd_russian_roulette_eps;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen_russian_roulette_eps(rd_russian_roulette_eps()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis_russian_roulette_eps(0 ,1);