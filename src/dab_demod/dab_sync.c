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
  - amplitude independent null symbol detection, SNR calculation
  - new phase reference symbol processing with combined time/frequency sync

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dab_sync.h"
#include "prs.h"

#define dbg 0

float snr2db( const float snr ) { return snr <= 1 ? 0 : 10 * log10(snr); }

void fft(    uint32_t size, fftwf_complex* in, fftwf_complex* out ) {
  fftwf_plan p = fftwf_plan_dft_1d( size, in, out, FFTW_FORWARD, FFTW_ESTIMATE );
  fftwf_execute(p);
  fftwf_destroy_plan(p);
}

void fftInv( uint32_t size, fftwf_complex* in, fftwf_complex* out ) {
  fftwf_plan p = fftwf_plan_dft_1d( size, in, out, FFTW_BACKWARD, FFTW_ESTIMATE );
  fftwf_execute(p);
  fftwf_destroy_plan(p);
}

uint32_t peakAmpIndex( uint32_t size, const fftwf_complex* vec ) {
  uint32_t maxPos = 0;
  float maxAmp = -1;
  for( uint32_t i = 0; i < size; ++i ) {
    const float amp = cabsf( vec[i] );
    if( amp > maxAmp ) {
      maxPos = i;
      maxAmp = amp;
    }
  }
  return maxPos;
}

/**
 * Corrects the frequency offset in [Hz] by multiplying with exp(-2Ipi*df*t*T).
 * The two TGUARD/2 before and after the TU data samples will also be corrected.
 * The phase reference symbol start before TGUARD/2 defines the time offset 0.
 */
void correctFrequency( dab_state* dab, float df, uint32_t size, const fftwf_complex* in, fftwf_complex* out ) {
  const float twoPiTdf = -2 * M_PI * df / SAMPLE_RATE;
  const fftwf_complex* prsStart = dab->frame + TNULL + dab->dt;
  const uint32_t indexOffset = in - prsStart; // Relative to PRS before GUARD/2.
  assert( indexOffset < TFRAME );
  const float phiOffset = twoPiTdf * indexOffset;
  const fftwf_complex corrStep = CMPLXF( cos(twoPiTdf), sin(twoPiTdf) );
  fftwf_complex corr = CMPLXF( cos(phiOffset), sin(phiOffset) );
  for( unsigned int i = 0; i < size; ++i ) {
    out[i] = in[i] * corr;
    corr *= corrStep;
  }
}

void gnuplot( const char* fn, uint32_t n, const fftwf_complex* ar ) {
#if dbg
  FILE* file = fopen( fn,"w+" );
  fprintf( file, "set border 3 # left + bottom\n" );
  fprintf( file, "set tics out\n" );
  fprintf( file, "set xtics nomirror\n" );
  fprintf( file, "set ytics nomirror\n" );
  fprintf( file, "set xzeroaxis\n" );
  fprintf( file, "set xlabel \"time [ms]\"\n" );
  fprintf( file, "set title \"%s\"\n", "title" );
  fprintf( file, "set style fill solid 0.2 noborder\n" );
  fprintf( file, "\n" );
  fprintf( file, "plot [0:%d] \\\n", n );
  fprintf( file, "  '-' using 1:2 notitle with lines lc 0 lt 1\n" );
  for( uint32_t k = 0; k < n; ++k )
    fprintf(file,"%.2f %.2f\n", (float)k, cabsf(ar[k]) );
  fprintf( file, "e\n" );
  fprintf( file, "\npause -1\n" );
  fclose(file);
#endif
}

/*
 * Find the null symbol start using a sliding window over the 1/16 sub-sampled
 * amplitude. Sets the null symbol start index * 2. Returns 0 if failed.
 * During the null symbol the signal is off, i.e. TNULL noise.
 * Compare this noise with TNULL/2 signal before and TNULL/2 signal after.
 */
int8_t dab_coarse_time_sync( dab_state* dab ) {
  const float WEAK_SIG = 4.0f;	// Weak signal amplitude.
  const float  LOW_SNR = 1.1f;	// Lower SNR threshold.
  const float HIGH_SNR = 2.0f;	// Upper SNR threshold.
  const uint32_t DEC = 16;	// Decimation = 1/16 sub-sampling.
  assert( TNULL % (2*DEC) == 0 );

  gnuplot( "frame.gp", TNULL+2*TS, dab->frame ); // First 3 symbols only

  // Fast coarse check of a matching prior time synchronization.
  if( dab->force_timesync == 0 ) {
    const float signal = ( cabsf(dab->frame[TNULL+  DEC])
			 + cabsf(dab->frame[TNULL+2*DEC])
			 + cabsf(dab->frame[TNULL+3*DEC])
			 + cabsf(dab->frame[TNULL+4*DEC]) ) / 4;
    const float  noise = ( cabsf(dab->frame[TNULL-  DEC])
			 + cabsf(dab->frame[TNULL-2*DEC])
			 + cabsf(dab->frame[TNULL-3*DEC])
			 + cabsf(dab->frame[TNULL-4*DEC]) ) / 4;
    const float snr = noise > 0 ? signal / noise : 1.0f;
    if( snr >= HIGH_SNR && signal >= WEAK_SIG ) {
      dab->dt = 0;
      dab->snr = snr;
      dab->signalAmp = signal;
      return -1; // Fasr check was successful.
    }
  }

  // Initialize sliding window with a null symbol stating at index 0.
  float amp[TFRAME/DEC];
  for( uint32_t i = 0; i < TFRAME/DEC; ++i ) amp[i] = cabsf(dab->frame[i*DEC]);
  float signalSum = 0;
  float noiseSum = 0;
  for( uint32_t i = 0; i < TNULL/DEC; ++i ) { // Null symbol = TNULL noise.
    noiseSum += amp[i];
  }
  for( uint32_t i = TNULL/DEC; i < 3*TNULL/(2*DEC); ++i ) { // Signal after.
    signalSum += amp[i];
  }
  const float noise = noiseSum / (TNULL/DEC);
  const float signal = signalSum / (TNULL/2/DEC);
  dab->snr = signal / noise;
  dab->signalAmp = signal;
  uint32_t bestIndex = 0;

  // Slide the window :   A--signal--B__noise__C--signal--D   -->
  for( uint32_t i = 3*TNULL/(2*DEC); i < TFRAME/DEC; ++i ) {
    const float ampD = amp[i];
    signalSum += ampD;
    const float ampC = amp[i-TNULL/(2*DEC)];
    signalSum -= ampC;
    noiseSum += ampC;
    const float ampB = amp[i-3*TNULL/(2*DEC)];
    noiseSum -= ampB;
    signalSum += ampB;
    uint32_t signalSamples;
    if( i*DEC >= 2*TNULL ) {
      const float ampA = amp[i-2*TNULL/DEC];
      signalSum -= ampA;
      signalSamples = TNULL/DEC;
    } else
      signalSamples = i-TNULL/DEC;

    // Calculate window S/N.
    const float signal = signalSum / signalSamples;
    const float noise =  noiseSum / (TNULL/DEC);
    const float snr =  signal / noise;
    if( snr > dab->snr ) { // Memorize best S/N and corresponding index.
      dab->snr = snr;
      dab->signalAmp = signal;
      bestIndex = i*DEC - TNULL - TNULL/2;
    }
  }

  // fprintf( stderr, "coarse_timeshift=%d\n", bestIndex );
  if( dab->snr < HIGH_SNR || dab->signalAmp < WEAK_SIG ) {
    dab->dt = 0;
    return 0; // No null symbol detected.
  } else if( bestIndex < TGUARD/2 && dab->force_timesync == 0 ) {
    dab->dt = 0; // Keep time synchronization as is.
    return -1;
  } else {
    dab->dt = bestIndex; // Null symbol with good S/N found.
    return -1;
  }

}

/**
 * Coarse frequency synchronization using the given TU carriers of the PRS FFT.
 * Should be done AFTER fine time and frequency synchronization.
 * Returns the PRS FFT spectrum shift in number of carriers (1kHz).
 * Zero is expected for a centered spectrum i.e. with channels [-768 .. 768].
 */
int32_t dab_coarse_freq_sync( const fftwf_complex* prsFFT ) {
  const int FREQ_HUB = 16; // +/-16kHz around the center frequency
  float bestAmp = 0;
  int32_t bestShift = 0;
  for( int k = -FREQ_HUB; k <= FREQ_HUB; ++k ) {
    float amp = cabsf(prsFFT[k-1+CARRIERS/2])
	      + cabsf(prsFFT[k  +CARRIERS/2])
      	      - cabsf(prsFFT[k+1+CARRIERS/2])
      	      - cabsf(prsFFT[k+2+CARRIERS/2])
	      - cabsf(prsFFT[k-2+TU-CARRIERS/2])
	      - cabsf(prsFFT[k-1+TU-CARRIERS/2])
	      + cabsf(prsFFT[k+  TU-CARRIERS/2])
	      + cabsf(prsFFT[k+1+TU-CARRIERS/2]);
    // fprintf( stderr, "fShift=%d x=%.1f\n", k, amp );
    if( amp > bestAmp ) {
      bestShift = k;
      bestAmp = amp;
    }
  }
  // fprintf( stderr, "bestShift=%d\n", bestShift );
  return bestShift;
}

/**
 * Fine frequency synchronization using the given TS samples in the time domain.
 * Calculates the phase difference of the GUARD/2 samples before and after the
 * TU symbol samples.
 */
float dab_fine_freq_sync( dab_state* dab, const fftwf_complex* symbol ) {
  const uint32_t NSYNC = TGUARD - 32;	// Use not the whole guard for sync.
  const fftwf_complex* p =		// Start of the guard used for sync.
    symbol + TGUARD - NSYNC;
  fftwf_complex mean = 0;
  for( uint32_t i = 0; i < NSYNC; ++i ) {
    mean += p[i] * conjf(p[i+TU]);
  }
  mean /= NSYNC;
  const float angle = cargf( mean );
  return -angle / ( 2 * M_PI ) * 1000;
}

/**
 * Phase reference symbol (PRS) correlation in the frequency domain.
 * Returns the synchronized PRS transformed to the frequency domain.
 * e.g. J.Cho "PC-based receiver for Eureka-147" 2001
 * e.g. K.Taura "A DAB receiver" 1996
 */
int8_t dab_fine_sync( dab_state* dab, fftwf_complex* prsFFT ) {
  fftwf_complex* prsCorr = fftwf_alloc_complex(TU); // Aligned malloc for SIMD.
  gnuplot( "frame.gp", TNULL+2*TS, dab->frame ); // First 3 symbols only

  // Correct frequency with dt estimated from prior frames.
  correctFrequency( dab, dab->df, TU, &dab->frame[TNULL+TGUARD/2+dab->dt], prsFFT );
  gnuplot( "prs.gp", TS, &dab->frame[TNULL] );

  //---------------------------------------------------------------------------
  // 1) Fine time synchronization by correlation of the PRS with the reference.
  //---------------------------------------------------------------------------
  fft( TU, prsFFT, prsFFT );
  int32_t dfCoarse = dab_coarse_freq_sync( prsFFT ); // To center FFT.
  gnuplot( "prs_fft.gp", TU, prsFFT );
  for( uint32_t k = 0; k < TU; ++k ) // Correlation = product of FFT's
    prsCorr[k] = prsFFT[(k+dfCoarse)%TU] * conjf(prs_ref[k]);
  gnuplot( "prs_corr_fft.gp", TU, prsCorr );
  fftInv( TU, prsCorr, prsCorr );
  gnuplot( "prs_corr.gp", TU, prsCorr );
  int32_t peakIndex = peakAmpIndex( TU, prsCorr );
  if( peakIndex >= TU/2 ) peakIndex -= TU;
  dab->dt = peakIndex - TGUARD/2; // Synchronize with ref before guard of PRS.
  // fprintf( stderr, "dt=%d\n", dab->dt );

  //---------------------------------------------------------------------------
  // 2) Fine frequency sync (must be done before coarse frequency sync).
  //---------------------------------------------------------------------------
  dab->df = dab_fine_freq_sync( dab, &dab->frame[TNULL+dab->dt] );
  correctFrequency( dab, dab->df, TU, &dab->frame[TNULL+TGUARD/2+dab->dt], prsFFT );

  //---------------------------------------------------------------------------
  // 3) Coarse frequency sync in units of one carrier = 1kHz.
  //---------------------------------------------------------------------------
  fft( TU, prsFFT, prsFFT );
  dfCoarse = dab_coarse_freq_sync( prsFFT );
  dab->df += dfCoarse * 1000.0f;
  if( dfCoarse ) { // Correct the PRS again using the corrected df.
    correctFrequency( dab, dab->df, TU, &dab->frame[TNULL+TGUARD/2+dab->dt], prsFFT );
    fft( TU, prsFFT, prsFFT );
  }

  // fprintf( stderr, "df=%d kHz\n", dfCoarse );

  //---------------------------------------------------------------------------
  // 4) Check PRS for correct time and frequency synchronization.
  //---------------------------------------------------------------------------
#if dbg
  fftwf_complex* prsTS = fftwf_alloc_complex(TS); // Aligned malloc for SIMD.
  correctFrequency( dab, dab->df, TS, &dab->frame[TNULL+dab->dt], prsTS );
  const float dfFine = dab_fine_freq_sync( dab, prsTS );
  fft( TU, prsTS + TGUARD/2, prsCorr );
  dfCoarse = dab_coarse_freq_sync( prsCorr );
  if( dfCoarse )
    fprintf( stderr, "Coarse frequency synchronization failed\n" );
  for( uint32_t k = 0; k < TU; ++k ) {
    prsCorr[k] *= conjf(prs_ref[k]);
  }
  gnuplot( "check_prs_corr_fft.gp", TU, prsCorr );
  fftInv( TU, prsCorr, prsCorr );
  gnuplot( "check_prs_corr.gp", TU, prsCorr );
  int32_t dt = peakAmpIndex( TU, prsCorr );
  if( dt >= TU/2 ) dt -= TU;
  dt -= TGUARD/2;
  const float df = dfFine + 1000.0f * dfCoarse;
  if( abs(dt) > 1 || fabsf(df) > 1.0f )
    fprintf( stderr, "sync check   dt=%d   df=%.1f\n", dt, df );
  fftwf_free(prsTS);
#endif

  fftwf_free(prsCorr);
  return -1; // success
}
