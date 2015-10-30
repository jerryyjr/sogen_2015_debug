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
structs.hpp contains every struct used in the program.
*/

#ifndef STRUCTS_HPP
#define STRUCTS_HPP

#include <cmath> // Needed for INFINITY
#include <cstdlib> // Needed for cmath
#include <cstring> // Needed for strlen, memset, memcpy
#include <iostream> // Needed for cout
#include <bitset> // Needed for bitset
#include <fstream> // Needed for ofstream

#include "macros.hpp"
#include "memory.hpp"

using namespace std;

char* copy_str(const char*); // init.h cannot be included because it requires this file, structs.h, creating a cyclical dependency; therefore, copy_str, declared in init.h, must be declared in this file as well in order to use it here

/* terminal contains colors, streams, and common messages for terminal output
	notes:
		There should be only one instance of terminal at any time.
	todo:
*/
struct terminal {
	// Escape codes
	const char* code_blue;
	const char* code_red;
	const char* code_reset;
	
	// Colors
	char* blue;
	char* red;
	char* reset;
	
	// Verbose stream
	streambuf* verbose_streambuf;
	ostream* verbose_stream;
	
	terminal () {
		this->code_blue = "\x1b[34m";
		this->code_red = "\x1b[31m";
		this->code_reset = "\x1b[0m";
		this->blue = copy_str(this->code_blue);
		this->red = copy_str(this->code_red);
		this->reset = copy_str(this->code_reset);
		this->verbose_stream = new ostream(cout.rdbuf());
	}
	
	~terminal () {
		mfree(this->blue);
		mfree(this->red);
		mfree(this->reset);
		delete verbose_stream;
	}
	
	// Indicates a task is done (pass terminal->verbose() into this function to print only with verbose mode on)
	void done (ostream& stream) {
		stream << this->blue << "Done" << this->reset << endl;
	}
	
	// Indicates a task is done
	void done () {
		done(cout);
	}
	
	// Indicates the program is out of memory
	void no_memory () {
		cout << this->red << "Not enough memory!" << this->reset << endl;
	}
	
	// Indicates the program couldn't read from a pipe
	void failed_pipe_read () {
		cout << this->red << "Couldn't read from the pipe!" << this->reset << endl;
	}
	
	// Indicates the program couldn't write to a pipe
	void failed_pipe_write () {
		cout << this->red << "Couldn't write to the pipe!" << this->reset << endl;
	}
	
	// Returns the verbose stream that prints only when verbose mode is on
	ostream& verbose () {
		return *(this->verbose_stream);
	}
	
	// Sets the stream buffer for verbose mode
	void set_verbose_streambuf (streambuf* sb) {
		this->verbose_stream->rdbuf(sb);
		this->verbose_streambuf = this->verbose_stream->rdbuf();
	}
};

/* input_params contains all of the program's input parameters (i.e. the given command-line arguments) as well as data associated with them
	notes:
		There should be only one instance of input_params at any time.
		Simulation state should not be kept in this struct; use sim_data instead.
		Variables should be initialized to the values indicated in the usage information.
	todo:
*/
struct input_params {
	// Input and output files' paths and names (either absolute or relative)
	char* params_file; // The path and name of the parameter sets file, default=input.params
	bool read_params; // Whether or not the read the parameter sets file, default=false
	char* ranges_file; // The path and name of the parameter ranges file, default=none
	bool read_ranges; // Whether or not to read the ranges file, default=false
	char* perturb_file; // The path and name of the perturbations file, default=none
	bool read_perturb; // Whether or not to read the perturbations file, default=false
	char* gradients_file; // The path and name of the gradients file, default=none
	bool read_gradients; // Whether or not to read the gradients file, default=false
	char* passed_file; // The path and name of the passed file, default=output.passed
	bool print_passed; // Whether or not to print the passed file, default=false
	char* dir_path; // The path of the output directory for concentrations or oscillation features, default=none
	bool print_cons; // Whether or not to print concentrations, default=false
	bool binary_cons_output; // Whether or not to print the binary or ASCII value of numbers in the concentrations output files
	char* features_file; // The path and file of the features file, default=none
	bool ant_features; // Whether or not to print oscillation features in the anterior
	bool post_features; // Whether or not to print oscillation features in the posterior
	bool print_features; // Whether or not to print general features file, default=false
	char* conditions_file; // The path and file of the conditions file, default=none
	bool print_conditions; // Whether or not to print which conditions passed, default=false
	char* scores_file; // The path and file of the scores file, default=none
	bool print_scores; // Whether or not to print the scores for every mutant, default=false
	int num_colls_print; // The number of columns of cells to print for plotting of single cells on top of each other
	
	// Sets
	int num_sets; // The number of parameter sets to simulate, default=1
	
	// Simulation size
	int width_total; // The width of the PSM, default=3
	int width_initial; // The width of the PSM before anterior growth, default=3
	int height; // The height of the tissue, default=1
	
	// Simulation details
	int time_total; // The number of minutes to run each simulation for, default=1200
	int time_split; // The number of minutes it takes for cells to split, default=6
	int time_til_growth; // The number of minutes to wait before allowing cells to grow into the anterior PSM, default=600
	int seed; // The seed, used for generating random numbers, default=generated from the time and process ID
	bool reset_seed; // Whether or not to reset the seed after each parameter set, default=false
	int pseed; // The seed, used for generating random parameter sets, default=generated from the time and process ID
	bool store_pseed; // Whether or not to store the parameter generation seed, pseed, default=false
	char* seed_file; // Default=none
	bool print_seeds; // Whether or not to print the seeds used to the seed file
	double step_size; // The time step in minutes used for Euler's method, default=0.01
	double max_con_thresh; // Maximum threshold for concentrations, default=INFINITY
	bool short_circuit; // Whether or not to stop simulating a parameter set after a mutant fails
	int num_active_mutants; // The number of mutants to simulate for each parameter set, default=num_mutants
	int big_gran; // The granularity in time steps with which to store data, default=1
	int small_gran; // The granularit in time steps with which to simulate data, default=1
	
	// Piping data
	bool piping; // Whether or not input and output should be piped (as opposed to written to disk), default=false
	int pipe_in; // The file descriptor to pipe data from, default=none (0)
	int pipe_out; // The file descriptor to pipe data into, default=none (0)
	
	// Output stream data
	bool verbose; // Whether or not the program is verbose, i.e. prints many messages about program and simulation state, default=false
	bool quiet; // Whether or not the program is quiet, i.e. redirects cout to /dev/null, default=false
	streambuf* cout_orig; // cout's original buffer to be restored at program completion
	ofstream* null_stream; // A stream to /dev/null that cout is redirected to if quiet mode is set
	
	input_params () {
		this->params_file = NULL;
		this->read_params = false;
		this->ranges_file = NULL;
		this->read_ranges = false;
		this->perturb_file = NULL;
		this->read_perturb = false;
		this->gradients_file = NULL;
		this->read_gradients = false;
		this->passed_file = NULL;
		this->print_passed = false;
		this->dir_path = NULL;
		this->print_cons = false;
		this->binary_cons_output = false;
		this->features_file = NULL;
		this->ant_features = false;
		this->post_features = false;
		this->print_features = false;
		this->conditions_file = NULL;
		this->print_conditions = false;
		this->scores_file = NULL;
		this->print_scores = false;
		this->num_colls_print = 0;
		this->num_sets = 1;
		this->big_gran = 1;
		this->small_gran = 1;
		this->width_total = 3;
		this->width_initial = 3;
		this->height = 1;
		this->time_total = 1200;
		this->time_split = 6;
		this->time_til_growth = 600;
		this->seed = 0;
		this->reset_seed = false;
		this->pseed = 0;
		this->store_pseed = false;
		this->seed_file = NULL;
		this->print_seeds = false;
		this->step_size = 0.01;
		this->max_con_thresh = INFINITY;
		this->short_circuit = false;
		this->num_active_mutants = NUM_MUTANTS;
		this->piping = false;
		this->pipe_in = 0;
		this->pipe_out = 0;
		this->verbose = false;
		this->quiet = false;
		this->cout_orig = NULL;
		this->null_stream = new ofstream("/dev/null");
	}
	
	~input_params () {
		mfree(this->params_file);
		mfree(this->perturb_file);
		mfree(this->gradients_file);
		mfree(this->passed_file);
		mfree(this->dir_path);
		mfree(this->features_file);
		mfree(this->conditions_file);
		mfree(this->scores_file);
		mfree(this->seed_file);
		delete this->null_stream;
	}
};

/* rates contains the rates specified by the current parameter set as well as perturbation and gradient data
	notes:
		There should be only one instance of rates at any time.
		rates_active is the final, active rates that should be used in the simulation.
	todo:
*/
struct rates {
	double rates_base[NUM_RATES]; // Base rates taken from the current parameter set
	double factors_perturb[NUM_RATES]; // Perturbations (as percentages with 1=100%) taken from the perturbations input file
	bool using_gradients; // Whether or not any rates have specified perturbations
	int width; // The total width of the simulation
	double* factors_gradient[NUM_RATES]; // Gradients (as arrays of (position, percentage with 1=100%) pairs) taken from the gradients input file
	bool has_gradient[NUM_RATES]; // Whether each rate has a specified gradient
	int cells; // The total number of cells in the simulation
	double* rates_cell[NUM_RATES]; // Rates per cell that factor in the base rates and each cell's perturbations
	double* rates_active[NUM_RATES]; // Rates per cell position that factor in the base rates, each cell's perburations, and the gradients at each position
	
	explicit rates (int width, int cells) {
		memset(this->rates_base, 0, sizeof(this->rates_base));
		memset(this->factors_perturb, 0, sizeof(this->factors_perturb));
		this->using_gradients = false;
		this->width = width;
		this->cells = cells;
		for (int i = 0; i < NUM_RATES; i++) {
			this->factors_gradient[i] = new double[width];
			for (int j = 0; j < width; j++) {
				this->factors_gradient[i][j] = 1;
			}
			this->has_gradient[i] = false;
			this->rates_cell[i] = new double[cells];
			this->rates_active[i] = new double[cells];
		}
	}
	
	~rates () {
		for (int i = 0; i < NUM_RATES; i++) {
			delete[] this->factors_gradient[i];
			delete[] this->rates_cell[i];
			delete[] this->rates_active[i];
		}
	}
};

/* con_levels contains concentration levels and active records for specific portions of a simulation
	notes:
		This is a general struct used in several places so make sure any changes are compatible with the main cl, baby_cl and each mutant's cl.
	todo:
*/
struct con_levels {
	bool initialized; // Whether or not this struct's data have been initialized
	int num_con_levels; // The number of concentration levels this struct stores (not necessarily the total number of concentration levels)
	int time_steps; // The number of time steps this struct stores concentrations for
	int cells; // The number of cells this struct stores concentrations for
	double*** cons; // A three dimensional array that stores [concentration levels][time steps][cells] in that order
	int* active_start_record; // Record of the start of the active PSM at each time step
	int* active_end_record; // Record of the end of the active PSM at each time step
	
	con_levels () {
		this->initialized = false;
	}
	
	con_levels (int num_con_levels, int time_steps, int cells, int active_start) {
		this->initialized = false;
		initialize(num_con_levels, time_steps, cells, active_start);
	}
	
	// Initializes the struct with the given number of concentration levels, time steps, and cells
	void initialize (int num_con_levels, int time_steps, int cells, int active_start) {
		// If the current size is big enough to fit the new size then reuse the memory, otherwise allocate the required memory
		if (this->initialized && this->num_con_levels >= num_con_levels && this->time_steps >= time_steps && this->cells >= cells) {
			this->reset();
			this->active_start_record[0] = active_start;
		} else {
			this->num_con_levels = num_con_levels;
			this->time_steps = time_steps;
			this->cells = cells;
			this->active_start_record = new int[time_steps];
			this->active_start_record[0] = active_start; // Initialize the active start record with the given position
			this->active_end_record = new int[time_steps];
			this->active_end_record[0] = 0; // Initialize the active end record at position 0
		
			this->cons = new double**[num_con_levels];
			for (int i = 0; i < num_con_levels; i++) {
				this->cons[i] = new double*[time_steps];
				for (int j = 0; j < time_steps; j++) {
					this->cons[i][j] = new double[cells];
					for (int k = 0; k < cells; k++) {
						this->cons[i][j][k] = 0; // Initialize every concentration level at every time step for every cell to 0
					}
				}
			}
			for (int j = 1; j < time_steps; j++) {
				this->active_start_record[j] = 0;
		        this->active_end_record[j] = 0;
			}
			this->initialized = true;
		}
	}
	
	// Sets every value in the struct to 0 but does not free any memory
	void reset () {
		if (this->initialized) {
			for (int i = 0; i < this->num_con_levels; i++) {
				for (int j = 0; j < this->time_steps; j++) {
					for (int k = 0; k < this->cells; k++) {
						this->cons[i][j][k] = 0;
					}
				}
			}
			for (int j = 0; j < this->time_steps; j++) {
				this->active_start_record[j] = 0;
				this->active_end_record[j] = 0;
			}
		}
	}
	
	// Frees the memory used by the struct
	void clear () {
		if (this->initialized) {
			for (int i = 0; i < this->num_con_levels; i++) {
				for (int j = 0; j < this->time_steps; j++) {
					delete[] this->cons[i][j];
				}
				delete[] this->cons[i];
			}
			delete[] this->cons;
			delete[] this->active_start_record;
            delete[] this->active_end_record;
			this->initialized = false;
		}
	}
	
	~con_levels () {
		this->clear();
	}
};

/* growin_array is an integer array that resizes when necessary by doubling its size until it can access the requested index
	notes:
		The [] operator has been overloaded for the sake of convenience. It returns a reference to the value requested but also resizes the array when necessary.
	todo:
*/
struct growin_array {
	int* array; // The array of integers
	int size; // The size of the array

	growin_array() {}

	void initialize(int size) {
		array = new int[size];
		this->size = size;
	}
	
	explicit growin_array (int size) {
		this->initialize(size);
	}

	int& operator[] (int index) {
		if (index >= this->size) {
			int* new_array = new int[2 * this->size];
			for (int i = 0; i < this->size; i++) {
				new_array[i] = array[i];
			}
			delete[] array;
			this->size = 2 * this->size;
			array = new_array;
		}
		return array[index];
	}
	
	void reset(int new_size) {
		this->size = new_size;
	}
	
	int get_size() {
		return this->size;
	}
	
	~growin_array () {
		delete[] array;
	}
};

/* features contains the oscillation features for a particular simulation
	notes:
	todo:
*/
struct features {
	double period_post[NUM_INDICES]; // The period of oscillations for relevant concentrations in the posterior
	double period_ant[NUM_INDICES]; // The period of oscillations for relevant concentrations in the anterior
	double amplitude_post[NUM_INDICES]; // The amplitude of oscillations for relevant concentrations in the posterior
	double amplitude_ant[NUM_INDICES]; // The amplitude of oscillations for relevant concentrations in the anterior
	double peaktotrough_mid[NUM_INDICES]; // The peak to trough ratio for relevant concentrations in the middle of the simulation time-wise
	double peaktotrough_end[NUM_INDICES]; // The peak to trough ratio for relevant concentrations at the end of the simulation time-wise
	double sync_score_post[NUM_INDICES]; // The synchronization score for the relevant concentrations in the posterior
	double sync_score_ant[NUM_INDICES]; // The synchronization score for the relevant concentrations in the anterior
	double num_good_somites[NUM_INDICES]; // The number of good somites for the relevant concentrations
	
	features () {
		memset(period_post, 0, sizeof(period_post));
		memset(period_ant, 0, sizeof(period_ant));
		memset(amplitude_post, 0, sizeof(amplitude_post));
		memset(amplitude_ant, 0, sizeof(amplitude_ant));
		memset(peaktotrough_mid, 0, sizeof(peaktotrough_mid));
		memset(peaktotrough_end, 0, sizeof(peaktotrough_end));
		memset(sync_score_post, 0, sizeof(sync_score_post));
		memset(num_good_somites, 0, sizeof(num_good_somites));
	}
};

/* mutant_data contains data for a particular mutant
	notes:
	todo:
*/
struct mutant_data {
	int index; // The index, i.e. how many mutants run before this one
	char* print_name; // The mutant's name for printing output
	char* dir_name; // The mutant's name for making its directory
	int num_knockouts; // The number of knockouts required to achieve this mutant
	int knockouts[2]; // The indices of the concentrations to knockout (num_knockouts determines how many indices to knockout)
	double overexpressions[NUM_INDICES]; // Overexpression factors with 1=100% overexpressed
	con_levels cl; // The concentration levels at the end of this mutant's posterior simulation run
	int (*tests[2])(mutant_data&, features&); // The posterior and anterior conditions tests
	int (*wave_test)(pair<int, int>[], int, mutant_data&, int, int); // The traveling wave conditions test
	int num_conditions[NUM_SECTIONS]; // The number of conditions this mutant is tested on
	int cond_scores[NUM_SECTIONS][MAX_CONDS_ANY]; // The score this mutant can achieve for each condition
	int max_cond_scores[NUM_SECTIONS]; // The maximum score this mutant can achieve for each section
	bool secs_passed[NUM_SECTIONS]; // Whether or not this mutant has passed each section's conditions
	bitset<1 + MAX_CONDS_ANY> conds_passed[NUM_SECTIONS]; // The score this mutant achieved for each condition when run
	features feat; // The oscillation features this mutant produced when run
	int print_con; // The index of the concentration that should be printed (usually mh1)
	
	mutant_data () {
		this->index = 0;
		this->print_name = NULL;
		this->dir_name = NULL;
		this->num_knockouts = 0;
		memset(this->knockouts, 0, sizeof(this->knockouts));
		memset(this->overexpressions, 0, sizeof(this->overexpressions));
		memset(this->tests, 0, sizeof(this->tests));
		this->wave_test = NULL;
		memset(this->num_conditions, 0, sizeof(this->num_conditions));
		memset(this->cond_scores, 0, sizeof(this->cond_scores));
		memset(this->max_cond_scores, 0, sizeof(this->max_cond_scores));
		memset(this->secs_passed, false, sizeof(this->secs_passed));
		this->print_con = CMH1;
	}
	
	~mutant_data () {
		mfree(this->print_name);
		mfree(this->dir_name);
		this->cl.clear();
	}
	
	// Calculates the maximum score this mutant can achieve for each section based on the scores given for each condition
	void calc_max_scores () {
		for (int i = 0; i < NUM_SECTIONS; i++) {
			this->max_cond_scores[i] = 0;
			for (int j = 0; j < this->num_conditions[i]; j++) {
				this->max_cond_scores[i] += this->cond_scores[i][j];
			}
		}
	}
};

/* sim_data contains simulation data, partially taken from input_params and partially derived from other information
	notes:
		There should be only one instance of sim_data at any time.
		sim_data copies some data from the input_params struct so simulation-irrelevant data does not have to be passed around in sim.cpp.
	todo:
*/
struct sim_data {
	// Times and timing
	double step_size; // The step size in minutes
	int time_total; // The number of minutes to run for
	int steps_total; // The number of time steps to simulate (total time / step size)
	int steps_split; // The number of time steps it takes for cells to split
	int steps_til_growth; // The number of time steps to wait before allowing cells to grow into the anterior PSM
	bool no_growth; // Whether or not the simulation should rerun with growth
	
	// Granularities
	int big_gran; // The granularity in time steps with which to analyze and store data
	int small_gran; // The granularit in time steps with which to simulate data
	
	// Cutoff values
	double max_con_thresh; // The maximum threshold concentrations can reach before the simulation is prematurely ended
	int max_delay_size; // The maximum number of time steps any delay in the current parameter set takes plus 1 (so that baby_cl and each mutant know how many minutes to store)
	
	// Sizes
	int width_total; // The width in cells of the PSM
	int width_initial; // The width in cells of the PSM before anterior growth
	int width_current; // The width in cells of the PSM at the current time step
	int height; // The height in cells of the PSM
	int cells_total; // The total number of cells of the PSM (total width * total height)
	
	// Neighbors and boundaries
	int** neighbors; // An array of neighbor indices for each cell position used in 2D simulations (2-cell and 1D calculate these on the fly)
	int active_start; // The start of the active portion of the PSM
	int active_end; // The end of the active portion of the PSM
	
	// PSM section and section-specific times
	int section; // Posterior or anterior (sec_post or sec_ant)
	int time_start; // The start time (in time steps) of the current simulation
	int time_end; // The end time (in time steps) of the current simulation
	int time_baby; // Time 0 for baby_cl at the end of a simulation
	
	// Mutants and condition scores
	int num_active_mutants; // The number of mutants to simulate for each parameter set
	int max_scores[NUM_SECTIONS]; // The maximum score possible for all mutants for each testing section
	int max_score_all; // The maximum score possible for all mutants for all testing sections
	
	explicit sim_data (input_params& ip) {
		this->step_size = ip.step_size;
		this->time_total = ip.time_total;
		this->steps_total = ip.time_total / ip.step_size;
		this->steps_split = ip.time_split / ip.step_size;
		this->steps_til_growth = ip.time_til_growth / ip.step_size;
		this->no_growth = this->steps_total == this->steps_til_growth || ip.width_initial == ip.width_total;
		this->big_gran = ip.big_gran;
		this->small_gran = ip.small_gran;
		this->max_con_thresh = ip.max_con_thresh;
		this->max_delay_size = 0;
		this->width_total = ip.width_total;
		this->width_initial = ip.width_initial;
		this->width_current = ip.width_initial;
		this->height = ip.height;
		this->cells_total = ip.width_total * ip.height;
		this->neighbors = new int*[this->cells_total];
		int num_neighbors;
		if (this->height == 1) {
			num_neighbors = NEIGHBORS_1D;
		} else {
			num_neighbors = NEIGHBORS_2D;
		}
		for (int k = 0; k < this->cells_total; k++) {
			this->neighbors[k] = new int[num_neighbors];
		}
		this->section = 0;
		this->time_start = 0;
		this->time_end = 0;
		this->time_baby = 0;
		this->num_active_mutants = ip.num_active_mutants;
		memset(this->max_scores, 0, sizeof(this->max_scores));
		this->max_score_all = 0;
	}
	
	// Initializes the scores once mutants have been initialized
	void initialize_conditions_data (mutant_data mds[]) {
		for (int i = 0; i < NUM_SECTIONS; i++) {
			for (int j = 0; j < this->num_active_mutants; j++) {
				this->max_scores[i] += mds[j].max_cond_scores[i];
			}
			this->max_score_all += this->max_scores[i];
		}
	}
	
	// Initializes the current width and the active positions before a simulation starts
	void initialize_active_data () {
		this->width_current = this->width_initial;
		this->active_start = this->width_initial - 1;
		this->active_end = (this->active_start - this->width_current + 1 + this->width_total) % this->width_total;
	}
	
	~sim_data () {
		for (int k = 0; k < this->cells_total; k++) {
			delete[] this->neighbors[k];
		}
		delete[] this->neighbors;
	}
};

/* input_data contains information for retrieving data from an input file
	notes:
		All input files should be read with read_file and an input_data struct, storing their contents in a string buffer.
	todo:
*/
struct input_data {
	char* filename; // The path and name of the file
	char* buffer; // A buffer to store the file's contents
	int size; // The number of bytes the file's contents take up
	int index; // The current index to access the buffer from
	
	explicit input_data (char* filename) {
		this->filename = filename;
		this->buffer = NULL;
		this->size = 0;
		this->index = 0;
	}
	
	~input_data () {
		mfree(this->buffer);
	}
};

/* st_context contains the spatiotemporal context at a particular point in the simulation
	notes:
	todo:
		TODO Rename this to something less horrid.
*/
struct st_context {
	int time_prev; // The previous time step
	int time_cur; // The current time step
	int cell; // The current cell
	
	explicit st_context (int time_prev, int time_cur, int cell) {
		this->time_prev = time_prev;
		this->time_cur = time_cur;
		this->cell = cell;
	}
};

/* di_args contains arguments for dim_int, the dimer interactions function
	notes:
		This struct is just a wrapper used to minimize the number of arguments passed into dim_int.
	todo:
*/
struct di_args {
	double** rs; // Active rates
	con_levels& cl; // Concentration levels
	st_context& stc; // Spatiotemporal context
	double* dimer_effects; // An array of dimer effects to store in which to store the results of dim_int
	
	explicit di_args (double** rs, con_levels& cl, st_context& stc, double dimer_effects[]) :
		rs(rs), cl(cl), stc(stc), dimer_effects(dimer_effects)
	{}
};

/* di_args contains indices for dim_int, the dimer interactions function
	notes:
		This struct is just a wrapper used to minimize the number of arguments passed into dim_int.
	todo:
*/
struct di_indices {
	int con_protein_self; // The index of this protein's concentrations
	int con_protein_other; // The index of the interacting protein's concentrations
	int con_dimer; // The index of this protein's and the interacting one's heterodimer's concentrations
	int rate_association; // The index of the heterodimer's rate of association
	int rate_dissociation; // The index of the heterodimer's rate of dissociation
	int dimer_effect; // The index in the dimer_effects array in the associated di_args struct
	
	explicit di_indices (int con_protein_self, int con_protein_other, int con_dimer, int rate_association, int rate_dissociation, int dimer_effect) :
		con_protein_self(con_protein_self), con_protein_other(con_protein_other), con_dimer(con_dimer), rate_association(rate_association), rate_dissociation(rate_dissociation), dimer_effect(dimer_effect)
	{}
};

/* cp_args contains indices for con_protein_her and con_protein_delta, the protein concentrations functions
	notes:
		This struct is just a wrapper used to minimize the number of arguments passed into con_protein_her and con_protein_delta.
	todo:
*/
struct cp_args {
	sim_data& sd; // Simulation data
	double** rs; // Active rates
	con_levels& cl; // Concentration levels
	st_context& stc; // Spatiotemporal context
	int* old_cells; // An array of cell indices at the start of each protein's delay
	double* dimer_effects; // An array of dimer effects calculated by dim_int
	
	explicit cp_args (sim_data& sd, double** rs, con_levels& cl, st_context& stc, int old_cells[], double dimer_effects[]) :
		sd(sd), rs(rs), cl(cl), stc(stc), old_cells(old_cells), dimer_effects(dimer_effects)
	{}
};

/* cph_indices contains indices for con_protein_her, the Her protein concentration function
	notes:
		This struct is just a wrapper used to minimize the number of arguments passed into con_protein_her.
	todo:
*/
struct cph_indices {
	int con_mrna; // The index of this protein's associated mRNA concentrations
	int con_protein; // The index of this protein's concentrations
	int con_dimer; // The index of this protein's homodimer's concentrations
	int rate_synthesis; // The index of this protein's rate of synthesis
	int rate_degradation; // The index of this protein's rate of degradation
	int rate_association; // The index of the homodimer's rate of association
	int rate_dissociation; // The index of the homodimer's rate of dissociation
	int delay_protein; // The index of this protein's rate of delay
	int dimer_effect; // The index in the dimer_effects array in the associated cp_args struct
	int old_cell; // The index in the old_cells array in the associated cp_args struct
	
	explicit cph_indices (int con_mrna, int con_protein, int con_dimer, int rate_synthesis, int rate_degradation, int rate_association, int rate_dissociation, int delay_protein, int dimer_effect, int old_cell) :
		con_mrna(con_mrna), con_protein(con_protein), con_dimer(con_dimer), rate_synthesis(rate_synthesis), rate_degradation(rate_degradation), rate_association(rate_association), rate_dissociation(rate_dissociation), delay_protein(delay_protein), dimer_effect(dimer_effect), old_cell(old_cell)
	{}
};

/* cpd_indices contains indices for con_protein_delta, the Delta protein concentration function
	notes:
		This struct is just a wrapper used to minimize the number of arguments passed into con_protein_delta.
	todo:
*/
struct cpd_indices {
	int con_mrna; // The index of this protein's associated mRNA concentrations
	int con_protein; // The index of this protein's concentrations
	int rate_synthesis; // The index of this protein's rate of synthesis
	int rate_degradation; // The index of this protein's rate of degradation
	int delay_protein; // The index of this protein's rate of delay
	int old_cell; // The index in the old_cells array in the associated cp_args struct
	
	explicit cpd_indices (int con_mrna, int con_protein, int rate_synthesis, int rate_degradation, int delay_protein, int old_cell) :
		con_mrna(con_mrna), con_protein(con_protein), rate_synthesis(rate_synthesis), rate_degradation(rate_degradation), delay_protein(delay_protein), old_cell(old_cell)
	{}
};

/* cd_args contains arguments for con_dimer, the dimer concentration function
	notes:
		This struct is just a wrapper used to minimize the number of arguments passed into con_dimer.
	todo:
*/
struct cd_args {
	sim_data& sd; // Simulation data
	double** rs; // Active rates
	con_levels& cl; // Concentration levels
	st_context& stc; // Spatiotemporal context
	
	explicit cd_args (sim_data& sd, double** rs, con_levels& cl, st_context& stc) :
		sd(sd), rs(rs), cl(cl), stc(stc)
	{}
};

/* cd_indices contains indices for con_dimer, the dimer concentration function
	notes:
		This struct is just a wrapper used to minimize the number of arguments passed into con_dimer.
	todo:
*/
struct cd_indices {
	int con_protein; // The index of this dimer's constituent protein concentration
	int rate_association; // The index of this dimer's rate of association
	int rate_dissociation; // The index of this dimer's rate of dissociation
	int rate_degradation; // The index of this dimer's rate of degradation
	
	explicit cd_indices (int con_protein, int rate_association, int rate_dissociation, int rate_degradation) :
		con_protein(con_protein), rate_association(rate_association), rate_dissociation(rate_dissociation), rate_degradation(rate_degradation)
	{}
};

#endif

