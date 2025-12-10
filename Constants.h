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
// Larger jump tables reduce average jump distance overhead and improve performance
// Current: 256 jumps = 24KB constant memory (3 arrays * 256 * 4 * 8 bytes)
#define NB_JUMP 256

// GPU group size (number of kangaroos processed together for batch modular inversion)
#define GPU_GRP_SIZE 128

// Optimal thread block size for modern GPUs (Ada/Blackwell)
// Used for __launch_bounds__ hint to compiler for register optimization
#define GPU_BLOCK_SIZE 256

// GPU number of runs per kernel call
// Higher values amortize kernel launch overhead and improve GPU occupancy
// Increased from 64 to 192 for better performance on modern GPUs
#define NB_RUN 192

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
