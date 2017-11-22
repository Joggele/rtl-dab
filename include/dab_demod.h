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
  - removed superfluous buffers because memcpy is expensive

*/

#pragma once

#include <complex.h>		// include before fftw3.h to use it there
#include <fftw3.h>
#include <stdint.h>
#include "dab_fifo.h"		// CircularBuffer

#define DEFAULT_BUF_LENGTH (16 * 16384)

typedef struct{

  // Buffers
  CircularBuffer fifo;		// Semaphore protected input buffer.
  fftwf_complex* frame;		// The complex I/Q samples of the whole frame.
  fftwf_complex* symbols[2];	// Current and previous symbol swap space.

  // Synchronization parameters.
  uint32_t frequency;		// Center carrier frequency in Hz.
  float df;			// Frequency shift in Hz.
  int32_t dt;			// Time shift in units of complex samples.
  int32_t startup_delay;
  uint8_t force_timesync;
  float signalAmp;		// Approximate amplitude measured in the PRS.
  float snr;			// Signal / noise ratio measured in null symbol.

  uint8_t symbols_demapped[75][1536*2]; // Demodulated and demapped symbols.
  uint8_t FIC[1536*2*3];	// FIC (3 symbols)
  uint8_t FIC_dep[3096*4];
  uint8_t FIC_dep_dec[768*4];
  uint8_t FIC_dep_dec_byte[(768*4)/8];

  uint8_t FIC_dep_enc[3096*4];  
  uint8_t FIC_enc_pun[1536*2*3];  

  uint8_t FIC_dep_dec_scr[768*4];
  uint8_t fib[12][256];
  //uint8_t fib_c[12][32];
  uint32_t f_interl_table[2048];

  /* fault injection */
  double p_e_prior_dep;
  double p_e_prior_vitdec;
  double p_e_after_vitdec;
  uint32_t bits_dab_frame_mantisse;

} dab_state;


int8_t dab_demod( dab_state* dab );
void dab_demod_init( dab_state* dab );
