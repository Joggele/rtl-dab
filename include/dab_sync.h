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
  - added DAB constants

*/

#include <math.h>
#include "dab_demod.h"

#define SAMPLE_RATE    	2048000	// Sample rate for TU = 2048 and 1kHz carriers.
#define TFRAME		 196608	// Transmission frame duration = 96ms.
#define TNULL		   2656	// Null symbol duration.
#define TS		   2552	// Symbol duration TU + GUARD.
#define TU		   2048	// 1ms symbol = 1kHz inverse carrier spacing.
#define TGUARD		    504	// Guard interval : TS - TU.
#define CARRIERS	   1536	// Number of transmitted carriers.

int8_t dab_coarse_time_sync( dab_state* dab );
int8_t dab_fine_sync( dab_state* dab, fftwf_complex* prsFFT );
int32_t dab_coarse_freq_sync( const fftwf_complex* prsFFT );
float dab_fine_freq_sync( dab_state* dab, const fftwf_complex* symbol );

void fft(    uint32_t size, fftwf_complex* in, fftwf_complex* out );
void fftInv( uint32_t size, fftwf_complex* in, fftwf_complex* out );
void correctFrequency( dab_state* dab, float df, uint32_t size, const fftwf_complex* in, fftwf_complex* out );
float snr2db( const float snr );
