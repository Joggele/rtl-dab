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
  - float fftwf is sufficient and much faster than double fftw
  - removed many superfluous buffers because memcpy is very expensive
  - new synchronization with software-based frequency correction
  - D-QPSK demodulation, deinterleaving and demapping together
  - dt breaks SIMD alignment inside a frame, so use aligned symbol for SIMD FFT

*/

#include <string.h>			// memcpy
#include "dab_demod.h"
#include "dab_sync.h"
#include "dab_fic_descramble.h"
#include "dab_fic_depuncture.h"
#include "dab_helper_functions.h"
#include "viterbi.h"

int8_t dab_demod( dab_state* dab ) {

  // Check for data in FIFO.
  if( dab->fifo.count < TFRAME*3 ) return 0;

  // Read input data from thread-safe FIFO.
  uint8_t buffer[2*TFRAME];
  if( dab->dt > 0 ) { // Skip dt samples and read a full frame.
    cbPop( &dab->fifo, 2*dab->dt );
    cbRead( &dab->fifo, buffer, 2*TFRAME );
  } else if( TFRAME + dab->dt > 0 ) { // Pad -dt samples and read rest frame.
    memset( buffer, 127, -2*dab->dt );
    cbRead( &dab->fifo, buffer - 2*dab->dt, 2*TFRAME + 2*dab->dt );
  }

  // Complex data conversion.
  for( uint32_t i = 0; i < TFRAME; ++i ) {
    dab->frame[i] = CMPLXF( buffer[i<<1] - 127, buffer[(i<<1)+1] - 127 );
  }

  // Coarse time synchronization using the null symbol.
  if( !dab_coarse_time_sync( dab ) ) { // No null symbol found or weak signal.
    fprintf( stderr, "snr =%5.1fdB   sigAmp =%3d   weak signal\n",
  		     snr2db(dab->snr), (int)dab->signalAmp );
    return 0;
  }
  dab->force_timesync = 0;
  if( dab->dt ) { // Not in sync so read remaining frame or go to next frame.
    fprintf( stderr, "snr =%5.1fdB   sigAmp =%3d   dt = %6d   resync\n",
  		     snr2db(dab->snr), (int)dab->signalAmp, dab->dt );
    if( dab->dt < 0 || dab->fifo.count < (uint32_t)dab->dt )
      return 0; // Skip this incomplete frame and wait for next.
    // Keep data after synchronization index and read the rest of this frame.
    memcpy( buffer, buffer + 2*dab->dt, 2*TFRAME - 2*dab->dt );
    cbRead( &dab->fifo, buffer + 2*TFRAME - 2*dab->dt, 2*dab->dt );
    dab->dt = 0; // Now the frame is complete and synchronized.
    for( uint32_t i = 0; i < TFRAME; ++i ) {
      dab->frame[i] = CMPLXF( buffer[i<<1] - 127, buffer[(i<<1)+1] - 127 );
    }
  }

  // Fine time and frequency synchronization using the phase reference symbol.
  fftwf_complex* previousSymbol = dab->symbols[0];
  fftwf_complex* currentSymbol = dab->symbols[1];
  if( !dab_fine_sync( dab, currentSymbol ) ) {
    fprintf( stderr, "snr =%5.1fdB   sigAmp =%3d   dt =%3d   sync failed\n",
		     snr2db(dab->snr), (int)dab->signalAmp, dab->dt );
    return 0;
  }
  fprintf( stderr, "snr =%5.1fdB   sigAmp =%3d   dt =%3d   df =%7.1fHz\n",
		   snr2db(dab->snr), (int)dab->signalAmp, dab->dt, dab->df );

  // Differential D-QPSK demodulation of the remaining 75 symbols.
  for( uint32_t i = 1; i < 76; ++i ) {

    // Swap previous and current symbol.
    fftwf_complex* tmp = previousSymbol;
    previousSymbol = currentSymbol;
    currentSymbol = tmp;

    // FFT transform symbol after frequency correction.
    correctFrequency( dab, dab->df, TU,
		      &dab->frame[TNULL+(TS*i)+TGUARD/2+dab->dt], currentSymbol );
    fft( TU, currentSymbol, currentSymbol );

    // D-QPSK symbol extraction.
    for( uint32_t j = 0; j < CARRIERS; ++j ) {

      // Frequency deinterleaving.
      const uint32_t k = dab->f_interl_table[j];

      // Differential D-QPSK demodulation : phase difference by complex division.
      const fftwf_complex sym = currentSymbol[k] * conjf(previousSymbol[k]);

      // Demapping
      dab->symbols_demapped[i][j]          = crealf(sym) > 0 ? 0 : 1;
      dab->symbols_demapped[i][CARRIERS+j] = cimagf(sym) > 0 ? 0 : 1;
    }

  }
  
  /* block partitioning */
  for (uint32_t i=0;i<CARRIERS*2;i++){
    dab->FIC[i+2*CARRIERS*0] = dab->symbols_demapped[1][i];
    dab->FIC[i+2*CARRIERS*1] = dab->symbols_demapped[2][i];
    dab->FIC[i+2*CARRIERS*2] = dab->symbols_demapped[3][i];
  }

  /* fault injection */
  if (dab->p_e_prior_dep>0) {
    binary_fault_injection(&dab->FIC[0],2304*4,dab->p_e_prior_dep);
  }
  
  /* FIC depuncture */
  dab_fic_depuncture(&dab->FIC[2304*0],&dab->FIC_dep[3096*0]);
  dab_fic_depuncture(&dab->FIC[2304*1],&dab->FIC_dep[3096*1]);
  dab_fic_depuncture(&dab->FIC[2304*2],&dab->FIC_dep[3096*2]);
  dab_fic_depuncture(&dab->FIC[2304*3],&dab->FIC_dep[3096*3]);
  
  /* fault injection */
  if (dab->p_e_prior_vitdec>0) {
    binary_fault_injection(&dab->FIC_dep[0],3096*4,dab->p_e_prior_vitdec);
  }

  /* Vitdec */
  viterbi( &dab->FIC_dep[3096*0], 3096, &dab->FIC_dep_dec[768*0]);
  viterbi( &dab->FIC_dep[3096*1], 3096, &dab->FIC_dep_dec[768*1]);
  viterbi( &dab->FIC_dep[3096*2], 3096, &dab->FIC_dep_dec[768*2]);
  viterbi( &dab->FIC_dep[3096*3], 3096, &dab->FIC_dep_dec[768*3]);

#if 0
  printf("fic\t");
  for (uint32_t i=0;i<3096;i++) printf( "%x", dab->FIC_dep[i] );
  printf("\n");
  printf("fic dec\t");
  for (uint32_t i=0;i<768;i++) printf( "%x", dab->FIC_dep_dec[i] );
  printf("\n");
#endif

  /* fault injection after vitdec 
     actually useless as frame has to be practically error free after viterbi 
     at least for DAB w/o Reed Solomon 
  */
  if (dab->p_e_after_vitdec>0) {
    binary_fault_injection(&dab->FIC_dep_dec[0],768*4,dab->p_e_after_vitdec);
  }

  /*De-scramble */
  dab_fic_descramble(&dab->FIC_dep_dec[768*0],&dab->FIC_dep_dec_scr[768*0],768);
  dab_fic_descramble(&dab->FIC_dep_dec[768*1],&dab->FIC_dep_dec_scr[768*1],768);
  dab_fic_descramble(&dab->FIC_dep_dec[768*2],&dab->FIC_dep_dec_scr[768*2],768);
  dab_fic_descramble(&dab->FIC_dep_dec[768*3],&dab->FIC_dep_dec_scr[768*3],768);

  /* FIC -> FIB */
  for (uint32_t i=0; i<256; i++) {
    dab->fib[0][i] = dab->FIC_dep_dec_scr[768*0+i];
    dab->fib[1][i] = dab->FIC_dep_dec_scr[768*0+256+i];
    dab->fib[2][i] = dab->FIC_dep_dec_scr[768*0+512+i];
    dab->fib[3][i] = dab->FIC_dep_dec_scr[768*1+i];
    dab->fib[4][i] = dab->FIC_dep_dec_scr[768*1+256+i];
    dab->fib[5][i] = dab->FIC_dep_dec_scr[768*1+512+i];
    dab->fib[6][i] = dab->FIC_dep_dec_scr[768*2+i];
    dab->fib[7][i] = dab->FIC_dep_dec_scr[768*2+256+i];
    dab->fib[8][i] = dab->FIC_dep_dec_scr[768*2+512+i];
    dab->fib[9][i] = dab->FIC_dep_dec_scr[768*3+i];
    dab->fib[10][i] = dab->FIC_dep_dec_scr[768*3+256+i];
    dab->fib[11][i] = dab->FIC_dep_dec_scr[768*3+512+i];
  }

return 1;

}

/* taken from openDAB */
void init_f_interl_table( dab_state* dab ) {
  int KI[TU];
  KI[0] = 0;
  for( int i = 1; i < TU; ++i ) {
    KI[i] = ( 13 * KI[i-1] + 511 ) % TU;
  }
  int n = 0;
  for( int i = 0; i < TU; ++i ) {
    if(  KI[i] >= 256 && KI[i] <= 1792 && KI[i] != 1024 ) {
      if( KI[i] < TU/2 ) {
	dab->f_interl_table[n++] = KI[i] + TU/2;
      } else {
        dab->f_interl_table[n++] = KI[i] - TU/2;
      }
    }
  }
}

void dab_demod_init(dab_state * dab){

  // Initialize synchronization parameters.
  dab->dt = 0;
  dab->df = 0;
  dab->startup_delay = 0;

  // Malloc of various buffers.
  cbInit( &dab->fifo, 4 * TFRAME*2 ); // circular buffer init for 4 frames
  dab->frame = fftwf_alloc_complex(TFRAME);
  dab->symbols[0] = fftwf_alloc_complex(TU); // Aligned malloc to use SIMD.
  dab->symbols[1] = fftwf_alloc_complex(TU); // Aligned malloc to use SIMD.
  init_f_interl_table(dab);	// init interleaver table
  init_viterbi();		// init viterbi decoder

  // Make sure to disable fault injection by default.
  dab->p_e_prior_dep = 0.0f;
  dab->p_e_prior_vitdec = 0.0f;
  dab->p_e_after_vitdec = 0.0f;
  dab->bits_dab_frame_mantisse = 52;
  
}
