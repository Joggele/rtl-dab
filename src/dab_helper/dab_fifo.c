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
  - repaired bottleneck using memcpy
  - semaphore protected to serve as input buffer
  - added cbPop

*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>	// memcpy
#include "dab_fifo.h"

/* http://en.wikipedia.org/wiki/Circular_buffer */

void cbInit( CircularBuffer* cb, uint32_t size ) {
  cb->size  = size;
  cb->start = 0;
  cb->count = 0;
  cb->elems = (uint8_t *)calloc(cb->size, sizeof(uint8_t));
  sem_init( &cb->ready, 0, 1 );
}

void cbFree( CircularBuffer* cb ) {
  sem_close( &cb->ready );
  free( cb->elems );
}

void lock( CircularBuffer* cb ) {
  sem_wait( &cb->ready ); // block if sem==0 else --sem
}

void unlock( CircularBuffer* cb ) {
  int val;
  sem_getvalue( &cb->ready, &val );
  if( !val ) sem_post( &cb->ready ); // ++sem
}

int cbIsFull( const CircularBuffer* cb) {
  return cb->count == cb->size;
}
 
int cbIsEmpty( const CircularBuffer* cb ) {
  return cb->count == 0;
}
 
void _cbPop( CircularBuffer* cb, uint32_t count ) {
  if( count > cb->count ) {
    fprintf( stderr, "sdr_fifo.cbPop : buffer underrun!\n" );
    cb->count = 0;
    cb->start = 0;
  } else {
    cb->count -= count;
    cb->start = ( cb->start + count ) % cb->size;
  }
}

void cbPop( CircularBuffer* cb, uint32_t count ) {
  lock(cb);
  _cbPop( cb, count );
  unlock(cb);
}

void cbWrite( CircularBuffer* cb, uint8_t* ar, uint32_t count ) {
  lock(cb);
  assert( cb->count <= cb->size );
  assert( cb->start <  cb->size );
  const uint32_t nFree = cb->size - cb->count;
  if( count > nFree ) {
    fprintf( stderr, "sdr_fifo.cbWrite : buffer overrun! skipping %d bytes\n",
	     count - nFree );
    _cbPop( cb, count - nFree );
  }
  uint32_t i = ( cb->start + cb->count ) % cb->size;
  const uint32_t nLeft = cb->size - i;
  if( count <= nLeft )
    memcpy( cb->elems + i, ar, count );
  else {
    memcpy( cb->elems + i, ar, nLeft );
    memcpy( cb->elems, ar + nLeft, count - nLeft );
  }
  cb->count += count;
  unlock(cb);
}
 
void cbRead( CircularBuffer* cb, uint8_t* ar, uint32_t count ) {
  lock(cb);
  assert( cb->count <= cb->size );
  assert( cb->start <  cb->size );
  if( count > cb->count ) {
    fprintf( stderr, "sdr_fifo.cbRead : buffer underrun! padding with %d zeros\n",
	     count - cb->count );
    memset( ar, 0, count - cb->count );
    ar += count - cb->count;
    count = cb->count;
  }
  const uint32_t nLeft = cb->size - cb->start;
  if( count <= nLeft )
    memcpy( ar, cb->elems + cb->start, count );
  else {
    memcpy( ar, cb->elems + cb->start, nLeft );
    memcpy( ar + nLeft, cb->elems, count - nLeft );
  }
  cb->count -= count;
  cb->start = ( cb->start + count ) % cb->size;
  unlock(cb);
}
