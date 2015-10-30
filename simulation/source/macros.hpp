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
macros.hpp contains every constant and macro function used in the program.
*/

#ifndef MACROS_HPP
#define MACROS_HPP

// Mutant indices
#define MUTANT_WILDTYPE		0
#define MUTANT_HER7			1
#define MUTANT_HER13		2
#define MUTANT_DELTA		3
#define MUTANT_HER7HER13	4
#define MUTANT_HER1			5
#define MUTANT_HER7DELTA	6
#define MUTANT_HER1DELTA	7
#define MUTANT_HER7OVER		8
#define MUTANT_HER1OVER		9
#define MUTANT_DELTAOVER	10

#define NUM_MUTANTS 		11

/// Named shortcuts for each concentration level of mRNA, protein, and dimer

// The time point at which the cell was born
#define BIRTH		0

// mRNA levels
#define CMH1		1
#define CMH7		2
#define CMH13		3
#define CMDELTA		4

// Protein levels
#define CPH1		5
#define CPH7		6
#define CPH13		7
#define CPDELTA		8

// Dimer levels
#define CPH1H1		9
#define CPH1H7		10
#define CPH1H13		11

#define CPH7H7		12
#define CPH7H13		13

#define CPH13H13	14

// The index of the cell's parent
#define PARENT		15

#define NUM_CON_LEVELS	29 // How big an array holding concentration levels must be
#define NUM_CON_STORE	7 // How many concentration levels to store in the bigger struct (for analyzing oscillation features)
#define MIN_CON_LEVEL	1 // The smallest index of a concetration level not including BIRTH or PARENT
#define MAX_CON_LEVEL	27 // The largest index of a concentration level not including BIRTH or PARENT

/// Named shortcuts for each rate of mRNA, protein, and dimer

// mRNA synthesis rates
#define RMSH1			0
#define RMSH7			1
#define RMSH13			2
#define RMSDELTA		3

// mRNA degradation rates
#define RMDH1			4
#define RMDH7			5
#define RMDH13			6
#define RMDDELTA		7

// Protein synthesis rates
#define RPSH1			8
#define RPSH7			9
#define RPSH13			10
#define RPSDELTA		11

// Protein degradation rates
#define RPDH1			12
#define RPDH7			13
#define RPDH13			14
#define RPDDELTA		15

// Dimer association rates
#define RDAH1H1			16
#define RDAH1H7			17
#define RDAH1H13		18
#define RDAH7H7			19
#define RDAH7H13		20
#define RDAH13H13		21

// Dimer dissociation rates
#define RDDIH1H1		22
#define RDDIH1H7		23
#define RDDIH1H13		24
#define RDDIH7H7		25
#define RDDIH7H13		26
#define RDDIH13H13		27

// Dimer degradation rates
#define RDDGH1H1		28
#define RDDGH1H7		29
#define RDDGH1H13		30
#define RDDGH7H7		31
#define RDDGH7H13		32
#define RDDGH13H13		33

// mRNA transcription delays
#define RDELAYMH1		34
#define RDELAYMH7		35
#define RDELAYMH13		36
#define RDELAYMDELTA	37

// Protein translation delays
#define RDELAYPH1		38
#define RDELAYPH7		39
#define RDELAYPH13		40
#define RDELAYPDELTA	41

// Critical number of molecules of proteins per cell for inhibition of transcription
#define RCRITPH1H1 		42
#define RCRITPH7H13		43

#define RCRITPDELTA		44

#define NUM_RATES		45 // How big an array holding rates must be
#define MIN_DELAY		34 // The smallest index referring to a delay
#define MAX_DELAY		41 // The largest index referring to a delay

/// Named shortcuts for each index of mRNA, protein, dimer

// Her gene indices
#define IH1		0
#define IH7		1
#define IH13	2

#define NUM_HER_INDICES	3 // How big an array holding the her gene indices must be

// mRNA indices
#define IMH1 	0
#define IMH7 	1
#define IMH13	2
#define IMDELTA	3

// Protein indices
#define IPH1	0
#define IPH7	1
#define IPH13	2
#define IPDELTA	3

#define NUM_INDICES		4 // How big an array holding mRNA or protein indices must be
#define NUM_DD_INDICES	2 // How many mRNA levels are affected by Delta and need to calculate its affects (stands for delta dependent indices)

// Neighbor counts
#define NEIGHBORS_1D	2
#define NEIGHBORS_2D	6

// Condition testing sections
#define SEC_POST		0
#define SEC_ANT			1
#define SEC_WAVE		2
#define NUM_SECTIONS	3

// Oscillation features
#define NUM_FEATS	2

// Condition weights
#define CW_A	5 // Weight for A class conditions
#define CW_B	0 // Weight for B class conditions
#define CW_AB	3 // Weight for A/B class conditions

// Maximum numbers of conditions
#define MAX_CONDS_POST	3
#define MAX_CONDS_ANT	3
#define MAX_CONDS_WAVE	4
#define MAX_CONDS_ALL	(MAX_CONDS_POST + MAX_CONDS_ANT + MAX_CONDS_WAVE)
#define MAX_CONDS_ANY	MAX(MAX(MAX_CONDS_POST, MAX_CONDS_ANT), MAX_CONDS_WAVE)

// Oscillation features
#define PERIOD			0
#define AMPLITUDE		1
#define SYNC			2
#define NUM_FEATURES	3
#define NUM_DATA_POINTS 10 // The number of data points required for synchronization plotting
#define INTERVAL 		60 // The length of the overlapping intervals for synchronization plotting

// Exit statuses
#define EXIT_SUCCESS			0
#define EXIT_MEMORY_ERROR		1
#define EXIT_FILE_READ_ERROR	2
#define EXIT_FILE_WRITE_ERROR	3
#define EXIT_PIPE_READ_ERROR	4
#define EXIT_PIPE_WRITE_ERROR	5
#define EXIT_INPUT_ERROR		6

// Macros for commonly used functions small enough to inject directly into the code
#define ABS(x) ((x) < 0 ? -(x) : (x))
#define SQUARE(x) ((x) * (x))
#define CUBE(x) ((x) * (x) * (x))
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) < (y) ? (y) : (x))
#define ZERO_IF_NAN(x) (is_nan(x) ? 0 : (x))
#define WRAP(x, y) ((x) + (y)) % (y)
#define INT_STRLEN(x) (log10((x) > 0 ? (x) : 1) + 1)

#endif

