/*	
	This file is part of the Snoopy code.

    Snoopy code is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Snoopy code is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Snoopy code.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "common.h"
#include "debug.h"
#include "libconfig/libconfig.h"

#define CONFIG_FILENAME		"snoopy.cfg"

read_config() {
	// Read the config file and initialize everyting
	config_t	config;		// Initialize the structure
	long tmp_v;
	
	const char * configname;
	
	DEBUG_START_FUNC;

	config_init(&config);
	
	if(!config_read_file(&config, CONFIG_FILENAME)) {
		MPI_Printf("Error reading configuration file in line %d: %s\n", config_error_line(&config), config_error_text(&config));
		ERROR_HANDLER(ERROR_CRITICAL, "Failed to read the configuration file");
	}

	if(config_lookup_string(&config,"configname",&configname)) {
		MPI_Printf("Using config file: %s.\n",configname);
	}
	// read physics parameters-------------------------------------------------------------------------------
	if(!config_lookup_float(&config, "physics.boxsize.[0]",&param.lx)) {
		param.lx = 1.0;
	}
	if(!config_lookup_float(&config, "physics.boxsize.[1]",&param.ly)) {
	param.ly = 1.0;
	}
	if(!config_lookup_float(&config, "physics.boxsize.[2]",&param.lz)) {
		param.lz = 1.0;
	}
	
	if(!config_lookup_float(&config, "physics.shear",&param.shear)) {
		param.shear = 0.0;
	}
	if(!config_lookup_float(&config, "physics.omega_shear",&param.omega_shear)) {
		param.omega_shear = 0.0;
	}
	if(!config_lookup_bool(&config, "code.antialiasing",&param.antialiasing)) {
		param.antialiasing = 1;
	}
	if(!config_lookup_float(&config, "output.timevar_step",&param.toutput_time)) {
		param.toutput_time = 1.0;
	}
	if(!config_lookup_float(&config, "output.snapshot_step",&param.toutput_flow)) {
		param.toutput_flow = 1.0;
	}

	config_destroy(&config);
	
	DEBUG_END_FUNC;
	
	return;
}
	
