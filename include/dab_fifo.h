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

#include <semaphore.h>
#include <stdint.h>

// JS TODO use C++ queue
typedef struct {
  uint32_t size;
  uint32_t start;
  uint32_t count;
  uint8_t* elems;
  sem_t    ready;
} CircularBuffer;

void	cbInit( CircularBuffer* cb, uint32_t size );
void	cbFree( CircularBuffer* cb );
int	cbIsFull( const CircularBuffer* cb );
int	cbIsEmpty( const CircularBuffer* cb );
void	cbPop( CircularBuffer* cb, uint32_t count );
void	cbWrite( CircularBuffer* cb, uint8_t* ar, uint32_t count );
void	cbRead( CircularBuffer* cb, uint8_t* ar, uint32_t count );
