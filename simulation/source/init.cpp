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
init.cpp contains initialization functions used before any simulations start.
Avoid adding biological functions in this file and add them to sim.cpp instead.
*/

#include <unistd.h> // Needed for getpid

#include "init.hpp" // Function declarations

#include "io.hpp"
#include "main.hpp"
#include "tests.hpp"

using namespace std;

terminal* term = NULL; // The global terminal struct

/* copy_str copies the given string, allocating enough memory for the new string
	parameters:
		str: the string to copy
	returns: a pointer to the new string
	notes:
	todo:
*/
char* copy_str (const char* str) {
	char* newstr = (char*)mallocate(sizeof(char) * strlen(str) + 1);
	return strcpy(newstr, str);
}

/* random_int generates a random integer in the range specififed by the given pair of integers
	parameters:
		range: a pair of integers that specifies the lower and upper bounds in that order
	returns: the random integer
	notes:
	todo:
*/
int random_int (pair<int, int> range) {
	return range.first + (rand() % (range.second - range.first + 1));
}

/* random_double generates a random double in the range specified by the given pair of doubles
	parameters:
		range: a pair of doubles that specifies the lower and upper bounds in that order
	returns: the random double
	notes:
	todo:
*/
double random_double (pair<double, double> range) {
	return range.first + (range.second - range.first) * rand() / (RAND_MAX + 1.0);
}

/* interpolate linearly interpolates the value at the given location between two given points
	parameters:
		x: the location at which to interpolate the value
		x0: the left location
		x1: the right location
		y0: the left location's value
		y1: the right location's value
	returns: the interpolated value
	notes:
		x0 and x1 are integers because this function is used to interpolate in a cell tissue.
	todo:
*/
inline double interpolate (double x, int x0, int x1, double y0, double y1) {
	return y0 + (y1 - y0) * ((x - x0) / (x1 - x0));
}

/* init_terminal creates and initializes a new terminal struct
	parameters:
	returns: nothing
	notes:
	todo:
*/
void init_terminal () {
	if (term != NULL) {
		delete term;
	}
	term = new terminal();
}

/* free_terminal frees the terminal from memory and resets the terminal text color to its default value
	parameters:
	returns: nothing
	notes:
	todo:
*/
void free_terminal () {
	cout << term->reset;
	delete term;
}

/* accept_input_params fills the given input_params with values from the given command-line arguments
	parameters:
		num_args: the number of command-line arguments (i.e. argc)
		args: the array of command-line arguments (i.e. argv)
		ip: the program's input parameters
	returns: nothing
	notes:
		Ensure the usage function's message in main.cpp matches the arguments this function accepts whenever editing or adding command-line argument acceptance.
	todo:
*/
void accept_input_params (int num_args, char** args, input_params& ip) {
	if (num_args > 1) { // If arguments were given (the 0th argument is the program name)
		for (int i = 1; i < num_args; i += 2) { // Iterate through each argument pair (if an argument does not have an accompanying value, i-- should be called when found)
			char* option = args[i];
			char* value;
			if (i < num_args - 1) {
				value = args[i + 1];
			} else {
				value = NULL;
			}
			
			// Accept command-line arguments in both short and long form
			if (option_set(option, "-i", "--params-file")) {
				ensure_nonempty(option, value);
				store_filename(&(ip.params_file), value);
				ip.read_params = true;
			} else if (option_set(option, "-R", "--ranges-file")) {
				ensure_nonempty(option, value);
				store_filename(&(ip.ranges_file), value);
				ip.read_ranges = true;
			} else if (option_set(option, "-A", "--anterior-feats")) {
				ip.ant_features = true;
				i--;				
			} else if (option_set(option, "-P", "--posterior-feats")) {
				ip.post_features = true;
				i--;
			} else if (option_set(option, "-u", "--perturb-file")) {
				ensure_nonempty(option, value);
				store_filename(&(ip.perturb_file), value);
				ip.read_perturb = true;
			} else if (option_set(option, "-r", "--gradients-file")) {
				ensure_nonempty(option, value);
				store_filename(&(ip.gradients_file), value);
				ip.read_gradients = true;
			} else if (option_set(option, "-o", "--print-passed")) {
				ensure_nonempty(option, value);
				store_filename(&(ip.passed_file), value);
				ip.print_passed = true;
			} else if (option_set(option, "-D", "--directory-path")) {
				ensure_nonempty(option, value);
				store_filename(&(ip.dir_path), value);
			} else if (option_set(option, "-t", "--print-cons")) {
				ip.print_cons = true;
				i--;
			} else if (option_set(option, "-B", "--binary-cons-output")) {
				ip.binary_cons_output = true;
				i--;
			} else if (option_set(option, "-f", "--print-osc-features")) {
				ensure_nonempty(option, value);
				store_filename(&(ip.features_file), value);
				ip.print_features = true;
			} else if (option_set(option, "-W", "--print-conditions")) {
				ensure_nonempty(option, value);
				store_filename(&(ip.conditions_file), value);
				ip.print_conditions = true;
			} else if (option_set(option, "-E", "--print-scores")) {
				ensure_nonempty(option, value);
				store_filename(&(ip.scores_file), value);
				ip.print_scores = true;
			} else if (option_set(option, "-L", "--print-cells")) {
				ensure_nonempty(option, value);
				ip.num_colls_print = atoi(value);
			} else if (option_set(option, "-b", "--big-granularity")) {
				ensure_nonempty(option, value);
				ip.big_gran = atoi(value);
				if (ip.big_gran < 1) {
					usage("The big granularity to simulate with must be a positive number of time steps. Set -b or --big-granularity to at least 1.");
				}
			} else if (option_set(option, "-g", "--small-granularity")) {
				ensure_nonempty(option, value);
				ip.small_gran = atoi(value);
				if (ip.small_gran < 1) {
					usage("The small granularity to simulate with must be a positive number of time steps. Set -g or --small-granularity to at least 1.");
				}
			} else if (option_set(option, "-x", "--total-width")) {
				ensure_nonempty(option, value);
				ip.width_total = atoi(value);
				if (ip.width_total < 2) {
					usage("The total width must be at least two cells. Set -x or --total-width to at least 2.");
				}
			} else if (option_set(option, "-w", "--initial-width")) {
				ensure_nonempty(option, value);
				ip.width_initial = atoi(value);
				if (ip.width_initial < 2) {
					usage("The initial width must be at least two cells. Set -w or --initial-width to at least 2.");
				}
			} else if (option_set(option, "-y", "--height")) {
				ensure_nonempty(option, value);
				ip.height = atoi(value);
				if (ip.height < 1) {
					usage("The tissue height must be at least one cell. Set -y or --height to at least 1.");
				}
			} else if (option_set(option, "-S", "--step-size")) {
				ensure_nonempty(option, value);
				ip.step_size = atof(value);
				if (ip.step_size <= 0) {
					usage("The time step for Euler's method must be a positive real number. Set -S or --time-step to be greater than 0.");
				}
			} else if (option_set(option, "-m", "--minutes")) {
				ensure_nonempty(option, value);
				ip.time_total = atoi(value);
				if (ip.time_total < 1) {
					usage("The simulation must be run for positive number of minutes. Set -m or --minutes to at least 1.");
				}
			} else if (option_set(option, "-T", "--split-time")) {
				ensure_nonempty(option, value);
				ip.time_split = atoi(value);
				if (ip.time_split < 0) {
					usage("Cells must split in a positive number of minutes (use 0 to stop cells from splitting). Set -T or --split-time 0 to at least 0.");
				} else if (ip.time_split == 0) {
					ip.time_split = (int)INFINITY;
				}
			} else if (option_set(option, "-G", "--time-til-growth")) {
				ensure_nonempty(option, value);
				ip.time_til_growth = atoi(value);
				if (ip.time_til_growth < 0) {
					usage("The time until cells are allowed to grow must be a nonnegative number of minutes. Set -G or --time-til-growth to at least 0.");
				}
			} else if (option_set(option, "-p", "--parameter-sets")) {
				ensure_nonempty(option, value);
				ip.num_sets = atoi(value);
				if (ip.num_sets < 1) {
					usage("The number of parameters to run must be a positive integer. Set -p or --parameter-sets to at least 1.");
				}
			} else if (option_set(option, "-s", "--seed")) {
				ensure_nonempty(option, value);
				ip.seed = atoi(value);
				if (ip.seed <= 0) {
					usage("The seed to generate random numbers must be a positive integer. Set -s or --seed to at least 1.");
				}
				ip.reset_seed = true;
			} else if (option_set(option, "-X", "--reset-seed")) {
				ip.reset_seed = true;
				i--;
			} else if (option_set(option, "-d", "--parameters-seed")) {
				ensure_nonempty(option, value);
				ip.pseed = atoi(value);
				ip.store_pseed = true;
				if (ip.pseed <= 0) {
					usage("The seed to generate random parameters must be a positive integer. Set -d or --parameters-seed to at least 1.");
				}
			} else if (option_set(option, "-e", "--print-seeds")) {
				ensure_nonempty(option, value);
				store_filename(&(ip.seed_file), value);
				ip.print_seeds = true;
			} else if (option_set(option, "-a", "--max-con-threshold")) {
				ensure_nonempty(option, value);
				ip.max_con_thresh = atof(value);
				if (ip.max_con_thresh < 0) {
					usage("The propensities threshold must be a nonnegative real number. Set -a or --max-con-threshold to be at least 0.");
				} else if (ip.max_con_thresh == 0) {
					ip.max_con_thresh = INFINITY;
				}
			} else if (option_set(option, "-C", "--short-circuit")) {
				ip.short_circuit = true;
				i--;
			} else if (option_set(option, "-M", "--mutants")) {
				ensure_nonempty(option, value);
				ip.num_active_mutants = atoi(value);
				if (ip.num_active_mutants < 1 || ip.num_active_mutants > NUM_MUTANTS) {
					int strlen_num = INT_STRLEN(NUM_MUTANTS);
					char* message = (char*)mallocate(sizeof(char) * (strlen("The number of mutants to run must be a positive integer up to the number of coded-in mutants. Set -M or --mutants to be at least 1 and no more than .") + strlen_num + 1));
					sprintf(message, "The number of mutants to run must be a positive integer up to the number of coded-in mutants. Set -M or --mutants to be at least 1 and no more than %d.", NUM_MUTANTS);
					usage(message);
				}
			} else if (option_set(option, "-I", "--pipe-in")) {
				ensure_nonempty(option, value);
				ip.piping = true;
				ip.pipe_in = atoi(value);
				if (ip.pipe_in <= 0) {
					usage("The file descriptor to pipe data from must be a positive integer. Set -I or --pipe-in to be at least 1.");
				}
			} else if (option_set(option, "-O", "--pipe-out")) {
				ensure_nonempty(option, value);
				ip.piping = true;
				ip.pipe_out = atoi(value);
				if (ip.pipe_out <= 0) {
					usage("The file descriptor to pipe data into must be a positive integer. Set -O or --pipe-out to be at least 1.");
				}
			} else if (option_set(option, "-c", "--no-color")) {
				mfree(term->blue);
				mfree(term->red);
				mfree(term->reset);
				strcpy(term->blue, "");
				strcpy(term->red, "");
				strcpy(term->reset, "");
				i--;
			} else if (option_set(option, "-v", "--verbose")) {
				if (!ip.verbose) {
					ip.verbose = true;
				}
				i--;
			} else if (option_set(option, "-q", "--quiet")) {
				if (!ip.quiet) {
					ip.quiet = true;
					ip.cout_orig = cout.rdbuf();
					cout.rdbuf(ip.null_stream->rdbuf());
					term->set_verbose_streambuf(ip.null_stream->rdbuf());
				}
				i--;
			} else if (option_set(option, "-h", "--help")) {
				usage("");
				i--;
			} else if (option_set(option, "-l", "--licensing")) { 
				licensing();
				i--;
			} else { // Do not ignore invalid arguments; exit with an error to let the user know this argument is problematic
				const char* message_0 = "'";
				const char* message_1 = "' is not a valid option! Please check that every argument matches one available in the following usage information.";
				char* message = (char*)mallocate(sizeof(char) * (strlen(message_0) + strlen(option) + strlen(message_1) + 1));
				sprintf(message, "%s%s%s", message_0, option, message_1);
				usage(message);
			}
		}
	}
}

/* option_set checks if the given string matches either given version (short or long) of an option
	parameters:
		option: the string to check
		short_name: the short version of the option
		long_name: the long version of the option
	returns: true if the string matches a version, false otherwise
	notes:
	todo:
*/
inline bool option_set (const char* option, const char* short_name, const char* long_name) {
	return strcmp(option, short_name) == 0 || strcmp(option, long_name) == 0;
}

/* ensure_nonempty ensures that an option that should have an associated value has one or exits with an error
	parameters:
		option: the option to check
		arg: the value to check for
	returns: nothing
	notes:
	todo:
*/
void ensure_nonempty (const char* option, const char* arg) {
	if (arg == NULL) {
		char* message = (char*)mallocate(strlen("Missing the argument for the '' option.") + strlen(option) + 1);
		sprintf(message, "Missing the argument for the '%s' option.", option);
		usage(message);
	}
}

/* check_input_params checks that the given command-line arguments are semantically valid
	parameters:
		ip: the program's input parameters
	returns: nothing
	notes:
		Since command-line arguments may be given in any order, it is impossible to check certain requirements immediately after an argument is read. This function is called after accept_input_params and therefore has access to every argument.
	todo:
*/
void check_input_params (input_params& ip) {
	if (ip.width_initial > ip.width_total) {
		usage("The initial width must be no more than the total width. Set the initial width (-w or --initial-width) to <= the total width (-x or --total-width).");
	}
	if (ip.height == 1 && ip.width_total < 2) {
		usage("The total width must be >= 2 for modeling cells. Set the total width (-x or --total-width) to >= 2.");
	}
	if (ip.height > 1 && (ip.width_total < 4 || ip.width_total % 2 != 0)) {
		usage("The total width must be >= 4 and even for modeling cell tissues (cell chains, i.e. height = 1, can have odd widths >= 3). Set the total width (-x or --total-width) to >= 4 and even or set the height (-y or --height) to 1.");
	}
	if (ip.big_gran < ip.small_gran) {
		usage("The big granularity must be at least the size of the small granularity. Set the big granularity (-b or --big-granularity) to >= the small granularity (-g or --small-granularity)");
	}
	if (ip.time_til_growth > ip.time_total) {
		usage("The time until growth must be no more than the total simulation time. Set the time until growth (-G or --time-til-growth) to <= the total time (-m or --total-time).");
	}
	if (ip.piping && (ip.pipe_in == 0 || ip.pipe_out == 0)) {
		usage("If one end of a pipe is specified, the other must be as well. Set the file descriptors for both the pipe in (-I or --pipe-in) and the pipe out (-O or --pipe-out).");
	}
	if (!(ip.piping || ip.read_params || ip.read_ranges)) {
		usage("Parameter must be piped in via -I or --pipe-in, read from a file via -i or --params-file, or generated from a ranges file and number of sets via -R or --ranges-file and -p or --parameter-sets, respectively.");
	}
	if (!ip.dir_path && (ip.print_cons || ip.ant_features || ip.post_features)) {
		usage("An output directory must be specified to print concentrations or oscillation features in the anterior.");
	}
	if ((ip.width_initial == ip.width_total || ip.time_til_growth == ip.time_total) && ip.ant_features) {
		usage("Printing oscillation features in the anterior was specified but there is not enough time for growth of the PSM in the anterior. Set the time until growth (-G or --time-til-growth) to < total time or set the initial width (-w or --initial-width) to less than the total width (-x or --total-width)");
	}
	if (!(ip.width_initial == ip.width_total || ip.time_til_growth == ip.time_total) && (ip.time_total < ip.time_til_growth + (ip.width_total - ip.width_initial) * ip.time_split + ip.width_total * ip.time_split)) {
		usage("Performing anterior simulations was specified but there is not enough time for the PSM to fill with cells at least twice. Set the total time (-m or --total-time) to longer.");
	}
	if (ip.reset_seed) {
		init_seeds(ip, 0, false, false);
	}
}

/* generate_seed generates a seed based on the current UNIX time and process ID (PID) of the program instance
	parameters:
	returns: the generated seed
	notes:
		Generating a seed based on both the time and PID ensures concurrently executed versions of this program do not generated the same seed.
	todo:
*/
int generate_seed () {
	return abs(((time(0) * 181) * ((getpid() - 83) * 359)) % 805306457);
}

/* init_seed initializes the simulation and parameter set seeds, printing each seed used to the seeds file if the user specified it
	parameters:
		ip: the program's input parameters
		set_num: the index of the set whose seed is being printed
		append: whether to append to the existing seeds file or create a new one, overwriting the old file
		indent_message: whether or not to indent the messages printed (based on when this function is called)
	returns: nothing
	notes:
	todo:
*/
void init_seeds (input_params& ip, int set_num, bool append, bool indent_message) {
	// If the file is being created then generate the parameter set seed
	if (!append && ip.store_pseed) {
		ip.pseed = generate_seed();
		if (indent_message) {
			term->verbose() << "  ";
		}
		term->verbose() << term->blue << "Using seed " << term->reset << ip.pseed << " for parameter set generation" << endl;
	}
	
	// If the seed was not specified by command-line argument or a new one should be generated, generate one
	if (ip.seed == 0 || append) {
		ip.seed = generate_seed();
		if (indent_message) {
			term->verbose() << "  ";
		}
		term->verbose() << term->blue << "Using seed " << term->reset << ip.seed << " for each run" << endl;
	}
	
	// If the user specified printing seeds to a file then print them to the specified file
	if (ip.print_seeds) {
		ofstream seed_file;
		term->verbose() << "  ";
		open_file(&seed_file, ip.seed_file, append);
		if (!append && ip.store_pseed) {
			seed_file << "pseed: " << ip.pseed << endl;
		}
		seed_file << "seed " << set_num << ": " << ip.seed << endl;
	}
}

/* reset_sesed resets the simulation seed and fast-forwards it if necessary (to avoid running an anterior simulation with the same random numbers as the previously run posterior)
	parameters:
		ip: the program's input parameters
		sd: the current simulation's data
	returns: nothing
	notes:
		If the simulation is changed to generated a different number of random numbers in the posterior then this function must be updated to reflect that.
	todo:
*/
void reset_seed (input_params& ip, sim_data& sd) {
	srand(ip.seed);
	if (sd.section == SEC_ANT) { // If simulating the anterior, return to the last location in the random number generation to avoid duplicate numbers
		for (int i = 0; i < NUM_RATES; i++) {
			for (int k = 0; k < sd.cells_total; k++) {
				rand();
			}
		}
	}
}

/* init_verbosity sets the verbose stream to /dev/null if verbose mode is not enabled
	parameters:
		ip: the program's input parameters
	returns: nothing
	notes:
	todo:
*/
void init_verbosity (input_params& ip) {
	if (!ip.verbose) {
		term->set_verbose_streambuf(ip.null_stream->rdbuf());
	}
}

/* read_sim_params fills in the parameter sets array via the method the user specified (piping from another program, a parameter sets file, or random generation from a ranges file)
	parameters:
		ip: the program's input parameters
		params_data: the input_data for the parameter sets input file
		sets: the array to add parameter sets to
		ranges_data: the input_data for the ranges input file
	returns: nothing
	notes:
		This function is responsible for filling in sets via whatever method the user specified so add any future input methods here.
	todo:
*/
void read_sim_params (input_params& ip, input_data& params_data, double**& sets, input_data& ranges_data) {
	cout << term->blue;
	if (ip.piping) { // If the user specified piping
		cout << "Reading pipe " << term->reset << "(file descriptor " << ip.pipe_in << ") . . . ";
		read_pipe(sets, ip);
		term->done();
	} else if (ip.read_params) { // If the user specified a parameter sets input file
		read_file(&params_data);
		sets = new double*[ip.num_sets];
		for (int i = 0; i < ip.num_sets; i++) {
			if (params_data.index < params_data.size) { // Parse only as many lines as specified, even if the file is longer
				sets[i] = new double[NUM_RATES];
				memset(sets[i], 0, sizeof(double) * NUM_RATES);
				if(!parse_param_line(sets[i], params_data.buffer, params_data.index)) { // Parse each line until the file is empty or the required number of sets have been found
					ip.num_sets = i;
				}
			} else {
				ip.num_sets = i;
			}
		}
	} else if (ip.read_ranges) { // If the user specified a ranges input file to generate random numbers from
		cout << "Generating " << term->reset << ip.num_sets << " random parameter sets according to the ranges in " << ranges_data.filename << " . . ." << endl;
		cout << "  ";
		read_file(&ranges_data);
		sets = new double*[ip.num_sets];
		pair <double, double> ranges[NUM_RATES];
		parse_ranges_file(ranges, ranges_data.buffer);
		srand(ip.pseed);
		for (int i = 0; i < ip.num_sets; i++) {
			sets[i] = new double[NUM_RATES];
			for (int j = 0; j < NUM_RATES; j++) {
				sets[i][j] = random_double(ranges[j]);
			}
		}
		term->done();
	} else {
		usage("Parameter must be piped in via -I or --pipe-in, read from a file via -i or --params-file, or generated from a ranges file and number of sets via -R or --ranges-file and -p or --parameter-sets, respectively.");
	}
}

/* read_perturb_params reads the data from the perturbations file if the user specified it
	parameters:
		ip: the program's input parameters
		perturb_data: the input_data for the perturbations input file
	returns: nothing
	notes:
	todo:
*/
void read_perturb_params (input_params& ip, input_data& perturb_data) {
	if (ip.read_perturb) {
		read_file(&perturb_data);
	}
}

/* read_gradients_params reads the data from the gradients file if the user specified it
	parameters:
		ip: the program's input parameters
		gradients_data: the input_data for the gradients input file
	returns: nothing
	notes:
	todo:
*/
void read_gradients_params (input_params& ip, input_data& gradients_data) {
	if (ip.read_gradients) {
		read_file(&gradients_data);
	}
}

/* fill_perturbations fills the factors_perturb array in the given rates struct based on the given perturbations input buffer
	parameters:
		rs: the current simulation's rates to fill
		perturbs: the perturbations buffer taken from the input file
	returns: nothing
	notes:
		The buffer should contain one perturbation factor per line in the format 'factor start end' where 'start' is the index of the concentration to start applying the perturbation and 'end' is the index of the concentration after which to stop applying the perturbation. A perturbation with a maximum absolute percentage of 'factor' is applied to every concentration from 'start' to 'end'.
	todo:
		TODO allow comments in files
*/
void fill_perturbations (rates& rs, char* perturbs) {
	if (perturbs != NULL) {
		static const char* usage_message = "There was an error reading the given perturbations file.";
		int index = 0;
		while (perturbs[index] != '\0') {
			double factor = 0; // The perturbation factor
			int con_start = 0; // The starting concentration
			int con_end = 0; // The ending concentration
			
			// Read the perturbations
			if (sscanf(perturbs + index, "%lf %d %d", &factor, &con_start, &con_end) == 3) {
				if (factor < 0) {
					usage("The given perturbations file includes at least one factor less than 0. Adjust the perturbations file given with -u or --perturb-file.");
				}
				if (con_start < 0 || con_start >= NUM_RATES || con_end < 0 || con_end >= NUM_RATES) {
					usage("The given perturbations file includes rates outside of the valid range. Adjust the perturbations file given with -u or --perturb-file or add the appropriate rates by editing the macros file and recompiling.");
				}
				
				// If the given factor and index range is valid then fill the current simulation's rates with the perturbations
				factor /= 100;
				for (int i = con_start; i <= con_end; i++) {
					rs.factors_perturb[i] = factor;
				}
			} else {
				usage(usage_message);
			}
			while (not_EOL(perturbs[index++])) {} // Skip to the next line
		}
	}
}

/* fill_gradients fills the factors_gradient array in the given rates struct based on the given gradients input buffer
	parameters:
		rs: the current simulation's rates to fill
		gradients: the gradients buffer taken from the input file
	returns: nothing
	notes:
		The buffer should contain one gradient per line in format 'concentration (position factor) (position factor) ...' with at least one (position factor) pair where 'concentration' is the index of the concentration to which to apply the gradient, 'position' is a column in the cell tissue, and 'factor' is the percentage the concentration should reach at the associated position.
	todo:
		TODO allow comments in files
*/
void fill_gradients (rates& rs, char* gradients) {
	if (gradients != NULL) {
		static const char* usage_message = "There was an error reading the given gradients file.";
		int con; // The index of the concentration
		int column; // The column in the cell tissue
		double factor; // The factor to apply
		int last_column = 0; // The last column in the cell tissue given a gradient factor
		int start_column; // The first column to apply the gradient from
		int i = 0; // The index in the buffer
		while (gradients[i] != '\0') {
			// Read the concentration value
			if (sscanf(gradients + i, "%d", &con) != 1) {
				usage(usage_message);
			}
			if (con < 0 || con > NUM_RATES) {
				usage("The given gradients file includes rate indices outside of the valid range. Please adjust the gradients file or add the appropriate rates by editing the macros file and recompiling.");
			}
			rs.using_gradients = true; // Mark that at least one concentration has a gradient
			rs.has_gradient[con] = true; // Mark that this concentration has a gradient
			
			// Read every (position factor) pair
			while (gradients[i++] != ' ') {} // Skip past the concentration index
			while (not_EOL(gradients[i])) {
				// Read a position factor pair
				if (sscanf(gradients + i, "(%d %lf)", &column, &factor) != 2) {
					usage(usage_message);
				}
				if (column < 0 || column >= rs.width) {
					usage("The given gradients file includes positions outside of the given simulation width. Please adjust the gradients file or increase the width of the simulation using -x or --total-width.");
				}
				if (factor < 0) {
					usage("The given gradients file includes factors less than 0. Please adjusted the gradients file.");
				}
				
				// Apply the gradient factor
				factor /= 100;
				start_column = last_column;
				last_column = column;
				for (int j = start_column + 1; j < column; j++) {
					rs.factors_gradient[con][j] = interpolate(j, start_column, column, rs.factors_gradient[con][start_column], factor);
				}
				rs.factors_gradient[con][column] = factor;
				while (gradients[i++] != ')') {} // Skip past the end of the pair
				while (gradients[i] == ' ') {i++;} // Skip any whitespace before the next pair
			}
			
			// Apply the last gradient factor to the rest of the columns
			for (int j = column + 1; j < rs.width; j++) {
				rs.factors_gradient[con][j] = rs.factors_gradient[con][column];
			}
			i++;
		}
	}
}

/* calc_max_delay_size calculates the maximum delay any given parameter set includes and returns that +1 to size con_levels structs
	parameters:
		ip: the program's input parameters
		sd: the current simulation's data
		rs: the current simulation's rates take perturbation factors from
		sets: the array of parameter sets to take delays from
	returns: nothing
	notes:
		This function calculates the maximum delay using every parameter set because this way con_levels structs that are sized based on the maximum delay do not have to be resized for every set.
	todo:
*/
void calc_max_delay_size (input_params& ip, sim_data& sd, rates& rs, double** sets) {
	double max = 0;
	for (int i = 0; i < ip.num_sets; i++) {
		for (int j = MIN_DELAY; j <= MAX_DELAY; j++) {
			for (int k = 0; k < sd.width_total; k++) {
				// Calculate the minimum delay, accounting for the maximum allowable perturbation and gradients
				max = MAX(max, (sets[i][j] + (sets[i][j] * rs.factors_perturb[j])) * rs.factors_gradient[j][k]);  
			}
		}
	}
	sd.max_delay_size = MIN(max, sd.time_total) / sd.step_size + 1; // If the maximum delay is longer than the simulation time then set the maximum delay to the simulation time
	if (sd.big_gran > sd.max_delay_size) {
		cout << term->red << "The given big granularity cannot be larger than the maximum delay time (in time steps) of any given parameter set! Please reduce the big granularity with -b or --big-granularity. Remember that adding perturbations to a delay will likely increase its duration." << term->reset << endl;
		exit(EXIT_INPUT_ERROR);
	}
}

/* delete_file closes the given file and frees it from memory
	parameters:
		file: a pointer to the output file stream to delete
	returns: nothing
	notes:
	todo:
*/
void delete_file (ofstream* file) {
	close_if_open(file);
	delete file;
}

/* create_passed_file creates a file to store parameter sets that passed
	parameters:
		ip: the program's input parameters
	returns: a pointer to the file output stream
	notes:
	todo:
*/
ofstream* create_passed_file (input_params& ip) {
	ofstream* file_passed = new ofstream();
	if (ip.print_passed) { // Print the passed sets only if the user specified it
		open_file(file_passed, ip.passed_file, false);
	}
	return file_passed;
}

/* create_cons_filenames creates the directories to store concentrations files for each mutant
	parameters:
		file: a pointer to the output file stream to delete
	returns: an array with each mutant's directory path
	notes:
	todo:
*/
char** create_dirs (input_params& ip, sim_data& sd, mutant_data mds[]) {
	char** dirnames = (char**)mallocate(sizeof(char*) * ip.num_active_mutants);
	if (ip.print_cons || ip.ant_features || ip.post_features) { // If printing concentrations or oscillation features in the anterior was specified by the user then create them
		// Find the path length
		int path_length = strlen(ip.dir_path);
		if (ip.dir_path[path_length - 1] == '/') { // Remove the trailing slash if necessary
			ip.dir_path[--path_length] = '\0';
		}
		
		// Create the output directory
		create_dir(ip.dir_path);
		
		// Allocate the required memory
		char* orig_dir_path = copy_str(ip.dir_path);
		size_t max_mut_length = 0;
		for (int i = 0; i < ip.num_active_mutants; i++) { // Calculate the length of the longest mutant string
			max_mut_length = MAX(max_mut_length, strlen(mds[i].dir_name));
		}
		char* full_dir_path = (char*)mallocate(sizeof(char) * (strlen(orig_dir_path) + max_mut_length + 3)); // +1 for each slash and 1 for the NULL terminator
		
		// Create the directories if required
		if (ip.print_cons || (ip.ant_features && !sd.no_growth) || ip.post_features) {
			// Create each mutant directory
			for (int i = 0; i < ip.num_active_mutants; i++) {
				strcpy(full_dir_path, orig_dir_path);
				strcat(full_dir_path, "/");
				strcat(full_dir_path, mds[i].dir_name);
				strcat(full_dir_path, "/");
				dirnames[i] = copy_str(full_dir_path);
				create_dir(dirnames[i]);
			}
		}
		
		// Free the allocated memory
		mfree(full_dir_path);
		mfree(orig_dir_path);
	} else { // If printing concentrations was not specified by the user then fill the array with NULL paths
		for (int i = 0; i < ip.num_active_mutants; i++) {
			dirnames[i] = NULL;
		}
	}
	return dirnames;
}

/* delete_cons_filenames frees the directory path strings from memory
	parameters:
		ip: the program's input parameters
		dirnames_cons: the array of mutant directory paths
	returns: nothing
	notes:
	todo:
*/
void delete_dirs (input_params& ip, char** dirnames) {
	if (dirnames != NULL) {
		for (int i = 0; i < ip.num_active_mutants; i++) {
			mfree(dirnames[i]);
		}
		mfree(dirnames);
	}
}

/* create_features_file creates a file to store the oscillation features of each simulation
	parameters:
		ip: the program's input parameters
		mds: the array of all mutant data
	returns: a pointer to the output file stream
	notes:
	todo:
*/
ofstream* create_features_file (input_params& ip, mutant_data mds[]) {
	ofstream* file_features = new ofstream();
	if (ip.print_features) { // Print the oscillation features only if the user specified it
		open_file(file_features, ip.features_file, false);
		
		// Print the file header
		*file_features << "set,";
		for (int i = 0; i < ip.num_active_mutants; i++) {
			*file_features << "post sync " << mds[i].print_name << ",post per " << mds[i].print_name << ",post amp " << mds[i].print_name << ",post per " <<  mds[i].print_name << "/wt,post amp " << mds[i].print_name << "/wt,";
			*file_features << "ant sync " << mds[i].print_name << ",ant per " << mds[i].print_name << ",ant amp " << mds[i].print_name << ",ant per " <<  mds[i].print_name << "/wt,ant amp " << mds[i].print_name << "/wt,";
		}
		*file_features << endl;
	}
	return file_features;
}

/* create_conditions_file creates a file to store the scores for each condition for each simulation
	parameters:
		ip: the program's input parameters
		mds: the array of all mutant data
	returns: a pointer to the output file stream
	notes:
	todo:
*/
ofstream* create_conditions_file (input_params& ip, mutant_data mds[]) {
	ofstream* file_conditions = new ofstream();
	if (ip.print_conditions) { // Print the condition scores only if the user specified it
		open_file(file_conditions, ip.conditions_file, false);
		
		// Print the file header
		*file_conditions << "set,";
		for (int i = 0; i < ip.num_active_mutants; i++) {
			*file_conditions << mds[i].print_name << ",";
			int index = 0;
			for (int j = 0; j < NUM_SECTIONS; j++) {
				for (int k = 0; k < mds[i].num_conditions[j]; k++) {
					*file_conditions << "cond" << index << ",";
					index++;
				}
			}
		}
		*file_conditions << "total score" << endl;
	}
	return file_conditions;
}

/* create_features_file creates a file to store the total scores for each mutant for each simulation
	parameters:
		ip: the program's input parameters
		mds: the array of all mutant data
	returns: a pointer to the output file stream
	notes:
	todo:
*/
ofstream* create_scores_file (input_params& ip, mutant_data mds[]) {
	ofstream* file_scores = new ofstream();
	if (ip.print_scores) { // Print the total scores only if the user specified it
		open_file(file_scores, ip.scores_file, false);
		
		// Print the file header
		*file_scores << "set,";
		for (int i = 0; i < ip.num_active_mutants; i++) {
			*file_scores << mds[i].print_name << " POST,WAVE,ANT,";
		}
		*file_scores << "Total Score" << endl;
	}
	return file_scores;
}

/* create_mutant_data fills in the name, knockouts, conditions, etc. for each mutant
	parameters:
		sd: the current simulation's data
	returns: the array of all mutant data
	notes:
	todo:
*/
mutant_data* create_mutant_data (sim_data& sd) {
	mutant_data* mds = new mutant_data[sd.num_active_mutants];
	
	// Index each mutant and initialize its concentration levels based on the maximum delay size
	for (int i = 0; i < sd.num_active_mutants; i++) {
		mds[i].index = i;
		mds[i].cl.initialize(NUM_CON_LEVELS, sd.max_delay_size, sd.cells_total, sd.active_start);
	}
	
	// Wild type
	mds[MUTANT_WILDTYPE].print_name = copy_str("wildtype");
	mds[MUTANT_WILDTYPE].dir_name = copy_str("wildtype");
	mds[MUTANT_WILDTYPE].num_knockouts = 0;
	mds[MUTANT_WILDTYPE].tests[SEC_POST] = test_wildtype_post;
	mds[MUTANT_WILDTYPE].tests[SEC_ANT] = test_wildtype_ant;
	mds[MUTANT_WILDTYPE].wave_test = test_wildtype_wave;
	mds[MUTANT_WILDTYPE].num_conditions[SEC_POST] = 3;
	mds[MUTANT_WILDTYPE].cond_scores[SEC_POST][0] = CW_A;
	mds[MUTANT_WILDTYPE].cond_scores[SEC_POST][1] = CW_A;
	mds[MUTANT_WILDTYPE].cond_scores[SEC_POST][2] = CW_A;
	mds[MUTANT_WILDTYPE].num_conditions[SEC_ANT] = 3;
	mds[MUTANT_WILDTYPE].cond_scores[SEC_ANT][0] = CW_A;
	mds[MUTANT_WILDTYPE].cond_scores[SEC_ANT][1] = CW_B;
	mds[MUTANT_WILDTYPE].cond_scores[SEC_ANT][2] = CW_B;
	mds[MUTANT_WILDTYPE].num_conditions[SEC_WAVE] = 4;
	mds[MUTANT_WILDTYPE].cond_scores[SEC_WAVE][0] = CW_B;
	mds[MUTANT_WILDTYPE].cond_scores[SEC_WAVE][1] = CW_B;
	mds[MUTANT_WILDTYPE].cond_scores[SEC_WAVE][2] = CW_B;
	mds[MUTANT_WILDTYPE].cond_scores[SEC_WAVE][3] = CW_B;
	mds[MUTANT_WILDTYPE].calc_max_scores();
	
	// Her7
	if (MUTANT_HER7 >= sd.num_active_mutants) {return mds;}
	mds[MUTANT_HER7].print_name = copy_str("her7 mutant");
	mds[MUTANT_HER7].dir_name = copy_str("her7");
	mds[MUTANT_HER7].num_knockouts = 1;
	mds[MUTANT_HER7].knockouts[0] = RPSH7;
	mds[MUTANT_HER7].tests[SEC_POST] = test_her7_mutant_post;
	mds[MUTANT_HER7].tests[SEC_ANT] = test_her7_mutant_ant;
	mds[MUTANT_HER7].num_conditions[SEC_POST] = 2;
	mds[MUTANT_HER7].cond_scores[SEC_POST][0] = CW_A;
	mds[MUTANT_HER7].cond_scores[SEC_POST][1] = CW_A;
	mds[MUTANT_HER7].num_conditions[SEC_ANT] = 2;
	mds[MUTANT_HER7].cond_scores[SEC_ANT][0] = CW_A;
	mds[MUTANT_HER7].cond_scores[SEC_ANT][1] = CW_B;
	mds[MUTANT_HER7].num_conditions[SEC_WAVE] = 0;
	mds[MUTANT_HER7].calc_max_scores();
	
	// Her13
	if (MUTANT_HER13 >= sd.num_active_mutants) {return mds;}
	mds[MUTANT_HER13].print_name = copy_str("her13 mutant");
	mds[MUTANT_HER13].dir_name = copy_str("her13");
	mds[MUTANT_HER13].num_knockouts = 1;
	mds[MUTANT_HER13].knockouts[0] = RPSH13;
	mds[MUTANT_HER13].tests[SEC_POST] = test_her13_mutant_post;
	mds[MUTANT_HER13].tests[SEC_ANT] = test_her13_mutant_ant;
	mds[MUTANT_HER13].num_conditions[SEC_POST] = 3;
	mds[MUTANT_HER13].cond_scores[SEC_POST][0] = CW_A;
	mds[MUTANT_HER13].cond_scores[SEC_POST][1] = CW_A;
	mds[MUTANT_HER13].cond_scores[SEC_POST][2] = CW_A;
	mds[MUTANT_HER13].num_conditions[SEC_ANT] = 2;
	mds[MUTANT_HER13].cond_scores[SEC_ANT][0] = CW_A;
	mds[MUTANT_HER13].cond_scores[SEC_ANT][1] = CW_B;
	mds[MUTANT_HER13].num_conditions[SEC_WAVE] = 0;
	mds[MUTANT_HER13].calc_max_scores();
	
	// Delta
	if (MUTANT_DELTA >= sd.num_active_mutants) {return mds;}
	mds[MUTANT_DELTA].print_name = copy_str("delta mutant");
	mds[MUTANT_DELTA].dir_name = copy_str("delta");
	mds[MUTANT_DELTA].num_knockouts = 1;
	mds[MUTANT_DELTA].knockouts[0] = RPSDELTA;
	mds[MUTANT_DELTA].tests[SEC_POST] = test_delta_mutant_post;
	mds[MUTANT_DELTA].tests[SEC_ANT] = test_delta_mutant_ant;
	mds[MUTANT_DELTA].num_conditions[SEC_POST] = 3;
	mds[MUTANT_DELTA].cond_scores[SEC_POST][0] = CW_A;
	mds[MUTANT_DELTA].cond_scores[SEC_POST][1] = CW_A;
	mds[MUTANT_DELTA].cond_scores[SEC_POST][2] = CW_A;
	mds[MUTANT_DELTA].num_conditions[SEC_ANT] = 2;
	mds[MUTANT_DELTA].cond_scores[SEC_ANT][0] = CW_A;
	mds[MUTANT_DELTA].cond_scores[SEC_ANT][1] = CW_A;
	mds[MUTANT_DELTA].num_conditions[SEC_WAVE] = 0;
	mds[MUTANT_DELTA].calc_max_scores();
	
	// Her7-Her13
	if (MUTANT_HER7HER13 >= sd.num_active_mutants) {return mds;}
	mds[MUTANT_HER7HER13].print_name = copy_str("her7-her13 mutant");
	mds[MUTANT_HER7HER13].dir_name = copy_str("her7her13");
	mds[MUTANT_HER7HER13].num_knockouts = 2;
	mds[MUTANT_HER7HER13].knockouts[0] = RPSH7;
	mds[MUTANT_HER7HER13].knockouts[1] = RPSH13;
	mds[MUTANT_HER7HER13].tests[SEC_POST] = test_her7her13_mutant_post;
	mds[MUTANT_HER7HER13].tests[SEC_ANT] = test_her7her13_mutant_ant;
	mds[MUTANT_HER7HER13].num_conditions[SEC_POST] = 3;
	mds[MUTANT_HER7HER13].cond_scores[SEC_POST][0] = CW_A;
	mds[MUTANT_HER7HER13].cond_scores[SEC_POST][1] = CW_A;
	mds[MUTANT_HER7HER13].cond_scores[SEC_POST][2] = CW_A;
	mds[MUTANT_HER7HER13].num_conditions[SEC_ANT] = 2;
	mds[MUTANT_HER7HER13].cond_scores[SEC_ANT][0] = CW_A;
	mds[MUTANT_HER7HER13].cond_scores[SEC_ANT][1] = CW_B;
	mds[MUTANT_HER7HER13].num_conditions[SEC_WAVE] = 0;
	mds[MUTANT_HER7HER13].calc_max_scores();
	
	// Her1
	if (MUTANT_HER1 >= sd.num_active_mutants) {return mds;}
	mds[MUTANT_HER1].print_name = copy_str("her1 mutant");
	mds[MUTANT_HER1].dir_name = copy_str("her1");
	mds[MUTANT_HER1].num_knockouts = 1;
	mds[MUTANT_HER1].knockouts[0] = RPSH1;
	mds[MUTANT_HER1].tests[SEC_POST] = test_her1_mutant_post;
	mds[MUTANT_HER1].tests[SEC_ANT] = test_her1_mutant_ant;
	mds[MUTANT_HER1].wave_test = test_her1_wave;
	mds[MUTANT_HER1].num_conditions[SEC_POST] = 3;
	mds[MUTANT_HER1].cond_scores[SEC_POST][0] = CW_A;
	mds[MUTANT_HER1].cond_scores[SEC_POST][1] = CW_A;
	mds[MUTANT_HER1].cond_scores[SEC_POST][2] = CW_A;
	mds[MUTANT_HER1].num_conditions[SEC_ANT] = 1;
	mds[MUTANT_HER1].cond_scores[SEC_ANT][0] = CW_A;
	mds[MUTANT_HER1].num_conditions[SEC_WAVE] = 1;
	mds[MUTANT_HER1].cond_scores[SEC_WAVE][0] = CW_AB;
	mds[MUTANT_HER1].calc_max_scores();
	
	// Her7-Delta
	if (MUTANT_HER7DELTA >= sd.num_active_mutants) {return mds;}
	mds[MUTANT_HER7DELTA].print_name = copy_str("her7-delta mutant");
	mds[MUTANT_HER7DELTA].dir_name = copy_str("her7delta");
	mds[MUTANT_HER7DELTA].num_knockouts = 2;
	mds[MUTANT_HER7DELTA].knockouts[0] = RPSH7;
	mds[MUTANT_HER7DELTA].knockouts[1] = RPSDELTA;
	mds[MUTANT_HER7DELTA].tests[SEC_POST] = test_her7delta_mutant_post;
	mds[MUTANT_HER7DELTA].tests[SEC_ANT] = test_her7delta_mutant_ant;
	mds[MUTANT_HER7DELTA].num_conditions[SEC_POST] = 1;
	mds[MUTANT_HER7DELTA].cond_scores[SEC_POST][0] = CW_AB;
	mds[MUTANT_HER7DELTA].num_conditions[SEC_ANT] = 2;
	mds[MUTANT_HER7DELTA].cond_scores[SEC_ANT][0] = CW_AB;
	mds[MUTANT_HER7DELTA].cond_scores[SEC_ANT][1] = CW_AB;
	mds[MUTANT_HER7DELTA].num_conditions[SEC_WAVE] = 0;
	mds[MUTANT_HER7DELTA].calc_max_scores();
	
	// Her1-Delta
	if (MUTANT_HER1DELTA >= sd.num_active_mutants) {return mds;}
	mds[MUTANT_HER1DELTA].print_name = copy_str("her1-delta mutant");
	mds[MUTANT_HER1DELTA].dir_name = copy_str("her7delta");
	mds[MUTANT_HER1DELTA].num_knockouts = 2;
	mds[MUTANT_HER1DELTA].knockouts[0] = RPSH1;
	mds[MUTANT_HER1DELTA].knockouts[1] = RPSDELTA;
	mds[MUTANT_HER1DELTA].tests[SEC_POST] = test_her1delta_mutant_post;
	mds[MUTANT_HER1DELTA].tests[SEC_ANT] = test_her1delta_mutant_ant;
	mds[MUTANT_HER1DELTA].num_conditions[SEC_POST] = 1;
	mds[MUTANT_HER1DELTA].cond_scores[SEC_POST][0] = CW_AB;
	mds[MUTANT_HER1DELTA].num_conditions[SEC_ANT] = 2;
	mds[MUTANT_HER1DELTA].cond_scores[SEC_ANT][0] = CW_AB;
	mds[MUTANT_HER1DELTA].cond_scores[SEC_ANT][1] = CW_AB;
	mds[MUTANT_HER1DELTA].num_conditions[SEC_WAVE] = 0;
	mds[MUTANT_HER1DELTA].calc_max_scores();
	
	// Her7-overexpressed
	if (MUTANT_HER7OVER >= sd.num_active_mutants) {return mds;}
	mds[MUTANT_HER7OVER].print_name = copy_str("her7-overexpressed mutant");
	mds[MUTANT_HER7OVER].dir_name = copy_str("her7over");
	mds[MUTANT_HER7OVER].num_knockouts = 0;
	mds[MUTANT_HER7OVER].overexpressions[RMSH7] = 2;
	mds[MUTANT_HER7OVER].tests[SEC_POST] = test_her7over_mutant_post;
	mds[MUTANT_HER7OVER].tests[SEC_ANT] = test_her7over_mutant_ant;
	mds[MUTANT_HER7OVER].num_conditions[SEC_POST] = 1;
	mds[MUTANT_HER7OVER].cond_scores[SEC_POST][0] = CW_A;
	mds[MUTANT_HER7OVER].num_conditions[SEC_ANT] = 1;
	mds[MUTANT_HER7OVER].cond_scores[SEC_ANT][0] = CW_AB;
	mds[MUTANT_HER7OVER].num_conditions[SEC_WAVE] = 0;
	mds[MUTANT_HER7OVER].calc_max_scores();
	
	// Her1-overexpressed
	if (MUTANT_HER1OVER >= sd.num_active_mutants) {return mds;}
	mds[MUTANT_HER1OVER].print_name = copy_str("her1-overexpressed mutant");
	mds[MUTANT_HER1OVER].dir_name = copy_str("her1over");
	mds[MUTANT_HER1OVER].num_knockouts = 0;
	mds[MUTANT_HER1OVER].overexpressions[RMSH1] = 2;
	mds[MUTANT_HER1OVER].tests[SEC_POST] = test_her1over_mutant_post;
	mds[MUTANT_HER1OVER].tests[SEC_ANT] = test_her1over_mutant_ant;
	mds[MUTANT_HER1OVER].num_conditions[SEC_POST] = 1;
	mds[MUTANT_HER1OVER].cond_scores[SEC_POST][0] = CW_A;
	mds[MUTANT_HER1OVER].num_conditions[SEC_ANT] = 1;
	mds[MUTANT_HER1OVER].cond_scores[SEC_ANT][0] = CW_AB;
	mds[MUTANT_HER1OVER].num_conditions[SEC_WAVE] = 0;
	mds[MUTANT_HER1OVER].calc_max_scores();
	mds[MUTANT_HER1OVER].print_con = CMH7;
	
	// Delta-overexpressed
	if (MUTANT_DELTAOVER >= sd.num_active_mutants) {return mds;}
	mds[MUTANT_DELTAOVER].print_name = copy_str("delta-overexpressed mutant");
	mds[MUTANT_DELTAOVER].dir_name = copy_str("deltaover");
	mds[MUTANT_DELTAOVER].num_knockouts = 0;
	mds[MUTANT_DELTAOVER].overexpressions[RMSDELTA] = 2;
	mds[MUTANT_DELTAOVER].tests[SEC_POST] = test_deltaover_mutant_post;
	mds[MUTANT_DELTAOVER].tests[SEC_ANT] = test_deltaover_mutant_ant;
	mds[MUTANT_DELTAOVER].num_conditions[SEC_POST] = 1;
	mds[MUTANT_DELTAOVER].cond_scores[SEC_POST][0] = CW_A;
	mds[MUTANT_DELTAOVER].num_conditions[SEC_ANT] = 2;
	mds[MUTANT_DELTAOVER].cond_scores[SEC_ANT][0] = CW_A;
	mds[MUTANT_DELTAOVER].cond_scores[SEC_ANT][1] = CW_AB;
	mds[MUTANT_DELTAOVER].num_conditions[SEC_WAVE] = 0;
	mds[MUTANT_DELTAOVER].calc_max_scores();
	
	return mds;
}

/* delete_mutant_data frees the given array of mutant data from memory
	parameters:
		mds: the array of mutant data
	returns: nothing
	notes:
	todo:
*/
void delete_mutant_data (mutant_data mds[]) {
	delete[] mds;
}

/* delete_sets frees the given array of parameter sets from memory
	parameters:
		sets: the array of parameter sets
		ip: the program's input parameters
	returns: nothing
	notes:
	todo:
*/
void delete_sets (double** sets, input_params& ip) {
	for (int i = 0; i < ip.num_sets; i++) {
		delete[] sets[i];
	}
	delete[] sets;
}

/* copy_cl_to_mutant copies the given concentration levels to the given mutant's concentration levels
	parameters:
		sd: the current simulation's data
		cl: the concentration levels to copy from
		md: the data of the mutant that just ran
	returns: nothing
	notes:
	todo:
*/
void copy_cl_to_mutant (sim_data& sd, con_levels& cl, mutant_data& md) {
	for (int i = 0; i < cl.num_con_levels; i++) {
		for (int j_sim = sd.time_baby, j_md = 0; j_md < sd.max_delay_size; j_sim = WRAP(j_sim + 1, sd.max_delay_size), j_md++) {
			for (int k = 0; k < cl.cells; k++) {
				md.cl.cons[i][j_md][k] = cl.cons[i][j_sim][k];
			}
			md.cl.active_start_record[j_md] = cl.active_start_record[j_sim];
			md.cl.active_end_record[j_md] = cl.active_end_record[j_sim];
		}
	}
}

/* copy_mutant_to_cl copies the concentration levels of the given mutant to the given concentration levels
	parameters:
		sd: the current simulation's data
		cl: the concentration levels to copy to
		md: the data of the mutant that will run again
	returns: nothing
	notes:
	todo:
*/
void copy_mutant_to_cl (sim_data& sd, con_levels& cl, mutant_data& md) {
	for (int i = 0; i < md.cl.num_con_levels; i++) {
		for (int j = 0; j < md.cl.time_steps; j++) {
			for (int k = 0; k < md.cl.cells; k++) {
				cl.cons[i][j][k] = md.cl.cons[i][j][k];
			}
			cl.active_start_record[j] = md.cl.active_start_record[j];
			cl.active_end_record[j] = md.cl.active_end_record[j];
		}
	}
	for (int j = 0; j < md.cl.time_steps; j++) {
		for (int k = 0; k < md.cl.cells; k++) {
			cl.cons[BIRTH][j][k] -= sd.steps_til_growth + sd.max_delay_size;
		}
	}
}

/* reset_cout resets the cout buffer to its original stream if quiet mode was on and cout was therefore redirected to /dev/null
	parameters:
		ip: the program's input parameters
	returns: nothing
	notes:
	todo:
*/
void reset_cout (input_params& ip) {
	if (ip.quiet) {
		cout.rdbuf(ip.cout_orig);
	}
}

