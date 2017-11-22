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

#include "dab_helper_functions.h"


void dab_bit_to_byte(uint8_t * in, uint8_t * out,uint32_t len) {
  uint32_t i;
  uint32_t bc=0;
  for (i=0;i<len;i+=8){
    out[bc] = (in[i+0]<<7) + (in[i+1]<<6) + (in[i+2]<<5) + (in[i+3]<<4) +   
      (in[i+4]<<3) + (in[i+5]<<2) + (in[i+6]<<1) +(in[i+7]<<0);
    bc++;
  }
}



int8_t dab_crc16(uint8_t * in,uint32_t len) {
  uint16_t POLY=0x8408;
  uint16_t POLY_CORRECT = 0xF0B8;
  uint16_t crc = 0xFFFF;
  uint16_t i;
  uint16_t c;
  for(i=0;i<len;i++){
    c = (crc & 0x0001) ^ ((in[i]>>0)&1);
    crc = crc >> 1;
    if (c)
      crc ^= POLY;
  }
  return (crc==POLY_CORRECT) ? 0 : 1;
}



void init_mt64(void) {
	uint64_t init[4]={UINT64_C(0x12345), UINT64_C(0x23456), UINT64_C(0x34567), UINT64_C(0x45678)}, length=4;
	init_by_array64(init, length);
}

uint8_t binary_fault_injection(uint8_t *in,uint32_t len,double p_e) {

  uint32_t i;
  double r;
  for (i=0;i<len;i++) {
    //r = ((double)rand())/RAND_MAX;
    r = genrand64_real2();
    in[i] = (r>p_e) ? in[i] : !in[i];
  }
  
  return 1;

}


uint8_t fftwf_complex_precision_reduction(fftwf_complex *in,uint32_t len,uint32_t precision) {
  
  uint32_t i,b;
  //  printf("%lf\n",in[0][0]);


  for (i=0;i<len;i++) {

    for (b=0;b<(52-precision);b++) {
      // unset bits
      ((uint32_t*)&in[i])[0] &= ~((uint64_t)1 << b);
      ((uint32_t*)&in[i])[1] &= ~((uint64_t)1 << b);
    }  
  }
  //printf("%lf\n",in[0][0]);

  return 1;
}

