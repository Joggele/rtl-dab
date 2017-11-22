/*
This file is part of rtl-dab
rtl-dab is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

rtl-dab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with rtl-dab.  If not, see <http://www.gnu.org/licenses/>.


david may 2012
david.may.muc@googlemail.com

*/

#include <errno.h>
#include <signal.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "dab_demod.h"
#include "dab_fic_parser.h"
#include "dab_analyzer.h"
#include "dab_helper_functions.h"



int corr_counter;

//ServiceInformation sinfo;
Ensemble sinfo;
Analyzer ana;

int main(int argc, char **argv){
  int i;
  dab_state dab;
  int frequency = 222055000;
  
  // open static file
  FILE *fh;

  // init mt
  init_mt64();

  // init demodulator structure
  dab_demod_init(&dab);

  // init FIC parser 
  dab_fic_parser_init(&sinfo);

  // init DAB Analyzer
  dab_analyzer_init(&ana);


  if (argc > 5) {
    fh = fopen(argv[1],"r");
    dab.bits_dab_frame_mantisse = atoi(argv[2]); 
    dab.p_e_prior_dep = atof(argv[3]);
    dab.p_e_prior_vitdec = atof(argv[4]);
    dab.p_e_after_vitdec = atof(argv[5]);
  } else {
    fh = fopen("/home/david/projects/rtl-dab/222055_dump.dump","r");
    dab.bits_dab_frame_mantisse = 52;
    dab.p_e_prior_dep = 0.0;
    dab.p_e_prior_vitdec = 0.0;
    dab.p_e_after_vitdec = 0.0;
  }


  // fill buffer and let the autogain settle
  uint8_t buffer[16*16384];
  for (i=0;i<20;i++) {
    unsigned int len = fread( buffer, 1, 16*16384, fh );
    cbWrite( &dab.fifo, buffer, len );
  }
  
  for (i=0;i<100;i++) {

    // read next dab frame
    unsigned int len = fread( buffer, 1, 16*16384, fh );
    cbWrite( &dab.fifo, buffer, len );

    // demodulate frame
    dab_demod(&dab);
    
    // parse FIC
    dab_fic_parser(dab.fib,&sinfo,&ana);
    
    // calculate error rates
    dab_analyzer_calculate_error_rates(&ana,&dab);
    
  }
  
  fprintf(stderr,"Analyzer: \n");
  fprintf(stderr,"received fibs: %i\n",ana.received_fibs);
  fprintf(stderr,"faulty   fibs: %i\n",ana.faulty_fibs);
  fprintf(stderr,"faulty fib rate: %f\n",(float)ana.faulty_fibs/(float)ana.received_fibs);
  fprintf(stderr,"mean channel ber: %f\n",ana.mean_ber);
  

  fclose(fh);
  
  fh = fopen("results.log","w+");
  fprintf(fh,"%f\n",ana.mean_ber);
  fprintf(fh,"%f",(float)ana.faulty_fibs/(float)ana.received_fibs);

  
  return 1;
}
