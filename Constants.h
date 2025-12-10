/*
* This file is part of the BSGS distribution (https://github.com/JeanLucPons/Kangaroo).
* Copyright (c) 2020 Jean Luc PONS.
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, version 3.
*
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
* General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CONSTANTSH
#define CONSTANTSH

// Release number
#define RELEASE "2.2"

// Use symmetry
//#define USE_SYMMETRY

// Number of random jumps
// Max 512 for the GPU (constant memory limit is 64KB)
// Increased conservatively to 64 for better jump distribution
#define NB_JUMP 64

// GPU group size (number of kangaroos processed together for batch modular inversion)
#define GPU_GRP_SIZE 128

// GPU number of runs per kernel call
// Increased conservatively from 64 to 96 for better kernel launch amortization
#define NB_RUN 96

// Kangaroo type
#define TAME 0  // Tame kangaroo
#define WILD 1  // Wild kangaroo

// SendDP Period in sec
#define SEND_PERIOD 2.0

// Timeout before closing connection idle client in sec
#define CLIENT_TIMEOUT 3600.0

// Number of merge partition
#define MERGE_PART 256

#endif //CONSTANTSH
