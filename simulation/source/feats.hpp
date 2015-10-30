/*
Simulation for zebrafish segmentation
Copyright (C) 2013 Ahmet Ay, Jack Holland, Adriana Sperlea, Sebastian Sangervasi

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
feats.hpp contains function declarations for feats.cpp.
*/

#ifndef FEATS_HPP
#define FEATS_HPP

#include "structs.hpp"

int get_peaks_and_troughs(sim_data&, con_levels&, int, int, growin_array&, growin_array&, growin_array&, int);
void osc_features_post(sim_data&, input_params&, con_levels&, features&, features&, char*, int, int, int, int);
void osc_features_ant(sim_data&, input_params&, features&, char*, con_levels&, mutant_data&, int, int, int, int, int);
int osc_features_wave(sim_data&, con_levels&, mutant_data&);
double post_sync(sim_data&, con_levels&, int, int, int);
double ant_sync(sim_data&, con_levels&, int, int);
void plot_ant_sync(sim_data&, con_levels&, int, ofstream*, bool);
double pearson_correlation(double*, double*, int, int);
int wave_testing(sim_data&, con_levels&, mutant_data&, int, int, int);

#endif
