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

Jörg Siegler 2017   dev dot js at web dot de
  - read raw samples from stdin, usage :
    'rtldab_test < rawSamples.dump' or 'cat rawSamples.dump | rtldab_test'

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

//ServiceInformation sinfo;
Ensemble sinfo;
Analyzer ana;
dab_state dab;

int main(void){
  int frequency = 222055000;
  
  // open static file
  FILE* fh = stdin;

  // init demodulator structure
  dab_demod_init(&dab);

  // init FIC parser 
  dab_fic_parser_init(&sinfo);

  // init DAB Analyzer
  dab_analyzer_init(&ana);

  uint8_t buffer[16*16384];
  while( 1 ) {

    // read next dab frame
    if( dab.fifo.count < 3*196608 ) {
      uint32_t len = fread( buffer, 1, 16*16384, fh );
      if( len <= 0 ) break;
      cbWrite( &dab.fifo, buffer, len );
    }

    // demodulate frame
    dab_demod(&dab);
    
    // parse FIC
    dab_fic_parser(dab.fib,&sinfo,&ana);

    // calculate error rates
    dab_analyzer_calculate_error_rates(&ana,&dab);

  } 
  
  fprintf(stderr,"ENSEMBLE STATUS:                                 \n");
  fprintf(stderr,"-------------------------------------------------\n");
  fprintf(stderr,"locked: %u \n",sinfo.locked);
  fprintf(stderr,"EnsembleLabel: %s \n",sinfo.esl->label);
  fprintf(stderr,"\nSubchannel Organization:                         \n");
  struct BasicSubchannelOrganization *sco;
  sco = sinfo.sco;
  while (sco->next != NULL) {
    fprintf(stderr,"SubChId: %2u | StartAddr: %4u | sl:%u | subchannelSize: %u   \n",
	    sco->SubChId,sco->startAddr,sco->shortlong,sco->subchannelSize);
    sco = sco->next;
  }
  fprintf(stderr,"\nService Information:                         \n");
  struct ServiceList *sl;
  sl = sinfo.sl;
  while (sl->next != NULL) {
    fprintf(stderr,"SId: %8X | SubChId: %2u | SCId %u\n",sl->SId,sl->scp->SubChId,sl->scp->SCId);
    sl = sl->next;
  }
  fprintf(stderr,"\nService Labels:                         \n");
  struct ProgrammeServiceLabel *psl;
  psl = sinfo.psl;
  while (psl->next != NULL) {
    fprintf(stderr,"SId: %8X | Label: %s \n",psl->SId,psl->label);
    psl = psl->next;
  }

  fprintf(stderr,"Analyzer: \n");
  fprintf(stderr,"received fibs: %i\n",ana.received_fibs);
  fprintf(stderr,"faulty   fibs: %i\n",ana.faulty_fibs);
  fprintf(stderr,"faulty fib rate: %f\n",(float)ana.faulty_fibs/(float)ana.received_fibs);
  fprintf(stderr,"mean channel ber: %f\n",ana.mean_ber);

  return 1;
}
