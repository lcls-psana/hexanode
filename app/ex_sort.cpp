#ifdef WIN32
	#pragma warning(disable : 4005)
	#pragma warning(disable : 4996)
#endif

#ifdef WIN64
	#pragma warning(disable : 4005)
	#pragma warning(disable : 4996)
#endif


#ifndef LINUX
	#ifndef WIN32
		#ifndef WIN64
			#define LINUX
		#endif
	#endif
#endif


#define NUM_IONS 200
#define NUM_CHANNELS 80

#include "math.h"
#include "hexanode_proxy/resort64c.h"
#include <stdio.h>
#include "stdlib.h"
#include "memory.h"



double	offset_sum_u, offset_sum_v, offset_sum_w;
double	w_offset, pos_offset_x, pos_offset_y;
int		command;
sort_class * sorter;





#ifndef LINUX
#include "conio.h"
__int32 my_kbhit()
{
	if (!_kbhit()) return 0;
	__int32 c = _getch();
	while (_kbhit()) c = _getch(); // empty keyboard buffer
	return c;
}
#else
#include <string.h>
#include <termios.h>
__int32 my_kbhit(void)
{
	struct termios term, oterm;
	__int32 fd = 0;
	__int32 c = 0;
	tcgetattr(fd, &oterm);
	memcpy(&term, &oterm, sizeof(term));
	term.c_lflag = term.c_lflag & (!ICANON);
	term.c_cc[VMIN] = 0;
	term.c_cc[VTIME] = 1;
	tcsetattr(fd, TCSANOW, &term);
	c = getchar();
	tcsetattr(fd, TCSANOW, &oterm);
	if (c == -1) return 0;
	return c;
}
#endif



void readline_from_config_file(FILE * ffile, char * text, __int32 max_len) {
	int i = -1;
	text[0] = 0;

	while (true) {
		i = -1;
		char c = 1;
		bool real_numbers_found = false;
		bool start_of_line = true;

		while (true) {
			fread(&c, 1, 1, ffile);
			if (c == 13) continue;
			if (start_of_line) {
				if (c==' ' || c==9 || c==10) continue;
				start_of_line = false;
			}
			if (c=='/') { // if there is a comment then read until end of the line
				while (c!=10) {
					fread(&c,1,1,ffile);
					if (c == 13) continue;
				}
				if (real_numbers_found) break;
				start_of_line = true;
				continue;
			}
			real_numbers_found = true;
			if (c==' ' || c==9) break;
			if (c == 10) break;
			if (c == ',') break;
			++i;
			if (i < max_len -1) {text[i] = c; text[i+1] = 0;}
		}
		if (real_numbers_found) break;
	}

	return;
}





int read_int(FILE * ffile) {
	char a[1024];
	readline_from_config_file(ffile, a, 1024);
	return atoi(a);
}



double read_double(FILE * ffile) {
	char a[1024];
	readline_from_config_file(ffile, a, 1024);
	return double(atof(a));
}



bool read_config_file(const char * name, sort_class *& sorter, int& command, double& offset_sum_u, double& offset_sum_v, double& offset_sum_w, double& w_offset, double& pos_offset_x, double& pos_offset_y)
{
	// read the config file:
	printf("opening %s... ",name);
	FILE * parameterfile_handle = fopen(name, "rt");

	if (!parameterfile_handle) {
		printf("file %s was not found.%c\n",name, 7);
		return false;
	}
	printf("ok\n");
	printf("reading %s... ",name);
	int int_dummy;
	command = read_int(parameterfile_handle);

	if (command == -1) {
		printf("ok\n");
		return false;
	}

	int_dummy = read_int(parameterfile_handle);
	sorter->use_HEX = int_dummy ? true: false;

	int_dummy = read_int(parameterfile_handle);
	sorter->common_start_mode = (!int_dummy) ? true:false;

	sorter->Cu1 = read_int(parameterfile_handle) -1;
	sorter->Cu2 = read_int(parameterfile_handle) -1;
	sorter->Cv1 = read_int(parameterfile_handle) -1;
	sorter->Cv2 = read_int(parameterfile_handle) -1;
	sorter->Cw1 = read_int(parameterfile_handle) -1;
	sorter->Cw2 = read_int(parameterfile_handle) -1;

	sorter->Cmcp   = read_int(parameterfile_handle) -1;
	sorter->use_MCP = (sorter->Cmcp > -1) ? true : false;

	offset_sum_u = read_double(parameterfile_handle);
	offset_sum_v = read_double(parameterfile_handle);
	offset_sum_w = read_double(parameterfile_handle);

	pos_offset_x = read_double(parameterfile_handle);
	pos_offset_y = read_double(parameterfile_handle);

	sorter->uncorrected_time_sum_half_width_u = read_double(parameterfile_handle);
	sorter->uncorrected_time_sum_half_width_v = read_double(parameterfile_handle);
	sorter->uncorrected_time_sum_half_width_w = read_double(parameterfile_handle);

	sorter->fu = 0.5*read_double(parameterfile_handle);
	sorter->fv = 0.5*read_double(parameterfile_handle);
	sorter->fw = 0.5*read_double(parameterfile_handle);
	w_offset   = read_double(parameterfile_handle);
	sorter->runtime_u = read_double(parameterfile_handle);
	sorter->runtime_v = read_double(parameterfile_handle);
	sorter->runtime_w = read_double(parameterfile_handle);
	
	sorter->MCP_radius  = read_double(parameterfile_handle);

	sorter->dead_time_anode = read_double(parameterfile_handle);
	sorter->dead_time_mcp   = read_double(parameterfile_handle);

	int_dummy = read_int(parameterfile_handle);
	sorter->use_sum_correction = (int_dummy != 0) ? true: false;
	int_dummy = read_int(parameterfile_handle);
	sorter->use_pos_correction = (int_dummy != 0) ? true: false;

	__int32 check_dummy = read_int(parameterfile_handle);
	if (check_dummy != 88888) {
		printf("File %s was not correctly read.\n", name);
		fclose(parameterfile_handle);
		if (sorter) {delete sorter; sorter = 0;}
		return false;
	}

	fclose(parameterfile_handle);
	// end of reading the config file:
	printf("ok\n");

	return true;
}







bool read_calibration_tables(const char * filename, sort_class * sorter)
{
	if (!filename) return false;
	if (!sorter) return false;

	FILE * infile_handle = fopen(filename,"rt");
	if (!infile_handle) return false;
	int points = 0;

	points = read_int(infile_handle);
	for (int j=0;j<points;++j) {
		double x = read_double(infile_handle);
		double y = read_double(infile_handle);
		if (sorter->use_sum_correction) sorter->signal_corrector->sum_corrector_U->set_point(x,y);
	}
	points = read_int(infile_handle);
	for (int j=0;j<points;++j) {
		double x = read_double(infile_handle);
		double y = read_double(infile_handle);
		if (sorter->use_sum_correction) sorter->signal_corrector->sum_corrector_V->set_point(x,y);
	}
	if (sorter->use_HEX) {
		points = read_int(infile_handle);
		for (int j=0;j<points;++j) {
			double x = read_double(infile_handle);
			double y = read_double(infile_handle);
			if (sorter->use_sum_correction) sorter->signal_corrector->sum_corrector_W->set_point(x,y);
		}
	}

	points = read_int(infile_handle);
	for (int j=0;j<points;++j) {
		double x = read_double(infile_handle);
		double y = read_double(infile_handle);
		if (sorter->use_pos_correction) sorter->signal_corrector->pos_corrector_U->set_point(x,y);
	}
	points = read_int(infile_handle);
	for (int j=0;j<points;++j) {
		double x = read_double(infile_handle);
		double y = read_double(infile_handle);
		if (sorter->use_pos_correction) sorter->signal_corrector->pos_corrector_V->set_point(x,y);
	}
	if (sorter->use_HEX) {
		points = read_int(infile_handle);
		for (int j=0;j<points;++j) {
			double x = read_double(infile_handle);
			double y = read_double(infile_handle);
			if (sorter->use_pos_correction) sorter->signal_corrector->pos_corrector_W->set_point(x,y);
		}
	}

	if (infile_handle) {fclose(infile_handle); infile_handle = 0;}
	return true;
}







bool create_calibration_tables(const char * filename, sort_class * sorter) 
{
	if (!sorter) return false;
	if (!filename) return false;
	FILE * fo = fopen(filename,"wt");
	sorter->do_calibration();
	int number_of_columns = sorter->sum_walk_calibrator->sumu_profile->number_of_columns;
	fprintf(fo,"\n\n%i  	// number of sum calibration points for layer U\n",number_of_columns);
	for (int binx=0; binx < number_of_columns; ++binx) {
		double x,y;
		sorter->sum_walk_calibrator->get_correction_point(x,y,binx,0); // 0 = layer u
		fprintf(fo,"%lg  %lg\n",x,y);
	}
	number_of_columns = sorter->sum_walk_calibrator->sumv_profile->number_of_columns;
	fprintf(fo,"\n\n%i  	// number of sum calibration points for layer V\n",number_of_columns);
	for (int binx=0; binx < number_of_columns; ++binx) {
		double x,y;
		sorter->sum_walk_calibrator->get_correction_point(x,y,binx,1); // 1 = layer v
		fprintf(fo,"%lg  %lg\n",x,y);
	}
	number_of_columns = sorter->sum_walk_calibrator->sumw_profile->number_of_columns;
	fprintf(fo,"\n\n%i  	// number of sum calibration points for layer W (only needed for HEX-detectors)\n",number_of_columns);
	if (sorter->use_HEX) {
		for (int binx=0; binx < number_of_columns; ++binx) {
			double x,y;
			sorter->sum_walk_calibrator->get_correction_point(x,y,binx,2); // 2 = layer w
			fprintf(fo,"%lg  %lg\n",x,y);
		}
	}

	number_of_columns = sorter->pos_walk_calibrator->number_of_columns;
	fprintf(fo,"\n\n%i  	// number of pos-calibration points for layer U\n",number_of_columns);
	for (int binx=0; binx < number_of_columns; ++binx) {
		double x,y;
		sorter->pos_walk_calibrator->get_correction_point(x,y,binx,0); // 0 = layer u
		fprintf(fo,"%lg  %lg\n",x,y);
	}
	number_of_columns = sorter->pos_walk_calibrator->number_of_columns;
	fprintf(fo,"\n\n%i  	// number of pos-calibration points for layer V\n",number_of_columns);
	for (int binx=0; binx < number_of_columns; ++binx) {
		double x,y;
		sorter->pos_walk_calibrator->get_correction_point(x,y,binx,1); // 1 = layer v
		fprintf(fo,"%lg  %lg\n",x,y);
	}
	number_of_columns = sorter->pos_walk_calibrator->number_of_columns;
	fprintf(fo,"\n\n%i  	// number of pos-calibration points for layer W (only needed for HEX-detectors)\n",number_of_columns);
	if (sorter->use_HEX) {
		for (int binx=0; binx < number_of_columns; ++binx) {
			double x,y;
			sorter->pos_walk_calibrator->get_correction_point(x,y,binx,2); // 2 = layer w
			fprintf(fo,"%lg  %lg\n",x,y);
		}
	}
	fclose(fo); fo = 0;
	return true;
}


















//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////







//////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[], char* envp[])
//////////////////////////////////////////////////////////////////////////////////////////
{
	printf("syntax: sort_LMF filename\n");
	printf("        This file will be sorted and\n");
	printf("        a new file will be written.\n\n");

	if (argc < 2)	{
		printf("Please provide a filename.\n");
		return 0;
	}
	if (argc > 2)	{
		printf("too many arguments\n");
		return 0;
	}

	double	tdc_ns[NUM_CHANNELS][NUM_IONS];
	__int32	number_of_hits[NUM_CHANNELS];
	command = -1;

	// The "command"-value is set in the first line of "sorter.txt"
	// 0 = only convert to new file format
	// 1 = sort and write new file 
	// 2 = calibrate fv, fw, w_offset
	// 3 = create calibration table files


	// create the sorter:
	sorter = new sort_class();

	if (!read_config_file("sorter.txt", sorter, command, offset_sum_u, offset_sum_v, offset_sum_w, w_offset, pos_offset_x, pos_offset_y)) {
		if (sorter) {delete sorter; sorter = 0;}
		return 0;
	}
	if (sorter) {
	if (sorter->use_sum_correction || sorter->use_pos_correction) {
		read_calibration_tables("calibration_table.txt", sorter);
	}
}
	
	if (command == -1) {
		printf("no config file was read. Nothing to do.\n");
		if (sorter) {delete sorter; sorter = 0;}
		return 0;
	}

	FILE * in_file = fopen(argv[1],"rb");

	if (!in_file) {
		printf("error: could not open the data file\n%s\n",argv[1]);
		return 0;
	}
	if (ferror(in_file)) {
		printf("error: could not open the data file\n%s\n",argv[1]);
		return 0;
	}

	// initialization of the sorter:
	printf("init sorter... ");
	sorter->TDC_resolution_ns = 0.025;
	sorter->tdc_array_row_length = NUM_IONS;
	sorter->count = (__int32*)number_of_hits;
	sorter->tdc_pointer = &tdc_ns[0][0];
	if (command >= 2) {
		sorter->create_scalefactors_calibrator(true, sorter->runtime_u, sorter->runtime_v, sorter->runtime_w, 0.78, sorter->fu, sorter->fv, sorter->fw); 
	}
 	int error_code = sorter->init_after_setting_parameters();
	if (error_code) {
		printf("sorter could not be initialized\n");
		char error_text[512];
		sorter->get_error_text(error_code,512,error_text);
		printf("Error %i: %s\n",error_code,error_text);
		return 0;
	}
	printf("ok\n");
	

	while (my_kbhit()); // empty keyboard buffer
	unsigned __int64 event_counter = 0;

	
	printf("reading event data... ");
	while (true) {
		event_counter++;
		if (event_counter%10000 == 0) {
			if (my_kbhit()) break;
		}

		//#error	TODO by end user:
		// Here you must read in a data block from your data file
		// and fill the array tdc_ns[][] and number_of_hits[]


		printf("error	Seaqrch for TODO by end user...");

		
		if (sorter->use_HEX) {
			// shift the time sums to zero:
			sorter->shift_sums(+1,offset_sum_u, offset_sum_v, offset_sum_w);
			// shift layer w so that the middle lines of all layers intersect in one point:
			sorter->shift_layer_w(+1,w_offset);
		} else {
			// shift the time sums to zero:
			sorter->shift_sums(+1,offset_sum_u, offset_sum_v);
		}
		// shift all signals from the anode so that the center of the detector is at x=y=0:
		sorter->shift_position_origin(+1, pos_offset_x, pos_offset_y);

		sorter->feed_calibration_data(true, w_offset); // for calibration of fv, fw, w_offset and correction tables
		if (sorter->scalefactors_calibrator) {
			if (sorter->scalefactors_calibrator->map_is_full_enough()) break;
		}
		
		int number_of_particles = 0;
		if (command == 1) {  // sort the TDC-Data and reconstruct missing signals
			// sorts/reconstructs the detector signals and apply the sum- and NL-correction.
			number_of_particles = sorter->sort();
			// "number_of_particles" is the number of reconstructed particles
		} else {
			number_of_particles = sorter->run_without_sorting();
		}
		


		printf("error	Seaqrch for TODO by end user...");
		//#error	TODO by end user:
		// write the results into a new data file.
		// the variable "number_of_particles" contains the number of
		// reconstructed particles.
		// the x and y  (in mm) and TOF (in ns) is stored in the array sorter->output_hit_array:
		// 
		// for the first particle:
		// sorter->output_hit_array[0]->x;
		// sorter->output_hit_array[0]->y;
		// sorter->output_hit_array[0]->time;
		// 
		// for the 2nd particle:
		// sorter->output_hit_array[1]->x;
		// sorter->output_hit_array[1]->y;
		// sorter->output_hit_array[1]->time;
		//
		// for each particle you can also retrieve the information about how the particle
		// was reconstructed (tog et some measure of the confidence):
		// sorter->output_hit_array[0]->method;


	} // end of the while loop


	if (command == 2) {
		printf("calibrating detector... ");
		sorter->do_calibration();
		printf("ok\n");
		if (sorter->scalefactors_calibrator) {
			printf("Good calibration factors are:\nf_U =%lg\nf_V =%lg\nf_W =%lg\nOffset on layer W =%lg\n", 2.*sorter->fu, 2.*sorter->scalefactors_calibrator->best_fv, 2.*sorter->scalefactors_calibrator->best_fw, sorter->scalefactors_calibrator->best_w_offset);
		}
	}


	if (command == 3) {   // generate and print correction tables for sum- and position-correction
		printf("creating calibration tables...\n");
		create_calibration_tables("calibration_table.txt", sorter);
		printf("\nfinished creating calibration tables\n");
	}

	
	if (sorter)	{delete sorter;  sorter  = 0;}
	

	return 0;
}
