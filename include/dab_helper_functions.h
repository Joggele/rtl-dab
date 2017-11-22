#include <stdint.h>
#include <stdlib.h>
#include <complex.h>	// include before fftw3.h to use it there
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>

#include <mt64.h>

void dab_bit_to_byte(uint8_t * in, uint8_t * out,uint32_t len);
int8_t dab_crc16(uint8_t * in,uint32_t len);

uint8_t binary_fault_injection(uint8_t *in,uint32_t len,double p_e);
uint8_t fftwf_complex_precision_reduction(fftwf_complex *in,uint32_t len,uint32_t precision);
void init_mt64(void);
