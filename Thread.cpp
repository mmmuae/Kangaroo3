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

#include "Kangaroo.h"
#include "Timer.h"
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>
#include <vector>
#ifndef WIN64
#include <pthread.h>
#endif

using namespace std;

// ----------------------------------------------------------------------------

#ifdef WIN64

THREAD_HANDLE Kangaroo::LaunchThread(LPTHREAD_START_ROUTINE func, TH_PARAM *p) {
  p->obj = this;
  return CreateThread(NULL, 0, func, (void*)(p), 0, NULL);
}
void  Kangaroo::JoinThreads(THREAD_HANDLE *handles,int nbThread) {
  WaitForMultipleObjects(nbThread, handles, TRUE, INFINITE);
}
void  Kangaroo::FreeHandles(THREAD_HANDLE *handles, int nbThread) {
  for (int i = 0; i < nbThread; i++)
    CloseHandle(handles[i]);
}
#else

THREAD_HANDLE Kangaroo::LaunchThread(void *(*func) (void *), TH_PARAM *p) {
  THREAD_HANDLE h;
  p->obj = this;
  pthread_create(&h, NULL, func, (void*)(p));
  return h;
}
void  Kangaroo::JoinThreads(THREAD_HANDLE *handles, int nbThread) {
  for (int i = 0; i < nbThread; i++)
    pthread_join(handles[i], NULL);
}
void  Kangaroo::FreeHandles(THREAD_HANDLE *handles, int nbThread) {
}
#endif

// ----------------------------------------------------------------------------

bool Kangaroo::isAlive(TH_PARAM *p) {

  bool isAlive = false;
  int total = nbCPUThread + nbGPUThread;
  for(int i=0;i<total;i++)
    isAlive = isAlive || p[i].isRunning;

  return isAlive;

}

// ----------------------------------------------------------------------------

bool Kangaroo::hasStarted(TH_PARAM *p) {

  bool hasStarted = true;
  int total = nbCPUThread + nbGPUThread;
  for (int i = 0; i < total; i++)
    hasStarted = hasStarted && p[i].hasStarted;

  return hasStarted;

}

// ----------------------------------------------------------------------------

bool Kangaroo::isWaiting(TH_PARAM *p) {

  bool isWaiting = true;
  int total = nbCPUThread + nbGPUThread;
  for (int i = 0; i < total; i++)
    isWaiting = isWaiting && p[i].isWaiting;

  return isWaiting;

}

// ----------------------------------------------------------------------------

uint64_t Kangaroo::getGPUCount() {

  uint64_t count = 0;
  for(int i = 0; i<nbGPUThread; i++)
    count += counters[0x80L + i];
  return count;

}

// ----------------------------------------------------------------------------

uint64_t Kangaroo::getCPUCount() {

  uint64_t count = 0;
  for(int i=0;i<nbCPUThread;i++)
    count += counters[i];
  return count;

}

// ----------------------------------------------------------------------------

string Kangaroo::GetTimeStr(double dTime) {

  char tmp[256];

  double nbDay = dTime / 86400.0;
  if (nbDay >= 1) {

    double nbYear = nbDay / 365.0;
    if (nbYear > 1) {
      if (nbYear < 5)
        sprintf(tmp, "%.1fy", nbYear);
      else
        sprintf(tmp, "%gy", nbYear);
    } else {
      sprintf(tmp, "%.1fd", nbDay);
    }

  } else {

    int iTime = (int)dTime;
    int nbHour = (int)((iTime % 86400) / 3600);
    int nbMin = (int)(((iTime % 86400) % 3600) / 60);
    int nbSec = (int)(iTime % 60);

    if (nbHour == 0) {
      if (nbMin == 0) {
        sprintf(tmp, "%02ds", nbSec);
      } else {
        sprintf(tmp, "%02d:%02d", nbMin, nbSec);
      }
    } else {
      sprintf(tmp, "%02d:%02d:%02d", nbHour, nbMin, nbSec);
    }

  }

  return string(tmp);

}

// ============================================================================
// LZB (Leading Zero Bits) Analysis Helper Functions
// ============================================================================

// Count leading zero bits in a 128-bit value
// Returns number of consecutive zero bits from MSB
uint32_t Kangaroo::CountLeadingZeroBits(int128_t *val) {

  uint32_t count = 0;

  // Check high 64 bits first
  if(val->i64[1] == 0) {
    count = 64;
    // High part is all zeros, check low part
    if(val->i64[0] == 0) {
      return 128; // Entire value is zero
    }
    // Count leading zeros in low part
    uint64_t v = val->i64[0];
    if(v == 0) return 128;

    // Binary search for leading zeros in 64-bit value
    if((v & 0xFFFFFFFF00000000ULL) == 0) { count += 32; v <<= 32; }
    if((v & 0xFFFF000000000000ULL) == 0) { count += 16; v <<= 16; }
    if((v & 0xFF00000000000000ULL) == 0) { count += 8;  v <<= 8;  }
    if((v & 0xF000000000000000ULL) == 0) { count += 4;  v <<= 4;  }
    if((v & 0xC000000000000000ULL) == 0) { count += 2;  v <<= 2;  }
    if((v & 0x8000000000000000ULL) == 0) { count += 1; }

  } else {
    // Count leading zeros in high part
    uint64_t v = val->i64[1];

    // Binary search for leading zeros in 64-bit value
    if((v & 0xFFFFFFFF00000000ULL) == 0) { count += 32; v <<= 32; }
    if((v & 0xFFFF000000000000ULL) == 0) { count += 16; v <<= 16; }
    if((v & 0xFF00000000000000ULL) == 0) { count += 8;  v <<= 8;  }
    if((v & 0xF000000000000000ULL) == 0) { count += 4;  v <<= 4;  }
    if((v & 0xC000000000000000ULL) == 0) { count += 2;  v <<= 2;  }
    if((v & 0x8000000000000000ULL) == 0) { count += 1; }
  }

  return count;
}

// Update the top 5 hottest buckets list
// Must be called with ghMutex locked
void Kangaroo::UpdateHotBuckets(uint32_t bucketId, uint32_t lzb, int128_t *xorDiff) {

  // Create new hot bucket entry
  HotBucket newBucket;
  newBucket.bucketId = bucketId;
  newBucket.lzb = lzb;
  newBucket.xorDiff.i64[0] = xorDiff->i64[0];
  newBucket.xorDiff.i64[1] = xorDiff->i64[1];

  // Check if this bucket is already in the list
  bool found = false;
  for(size_t i = 0; i < topHotBuckets.size(); i++) {
    if(topHotBuckets[i].bucketId == bucketId) {
      // Update if this is a better (higher) LZB for this bucket
      if(lzb > topHotBuckets[i].lzb) {
        topHotBuckets[i].lzb = lzb;
        topHotBuckets[i].xorDiff.i64[0] = xorDiff->i64[0];
        topHotBuckets[i].xorDiff.i64[1] = xorDiff->i64[1];
      }
      found = true;
      break;
    }
  }

  if(!found) {
    // Add new bucket
    topHotBuckets.push_back(newBucket);
  }

  // Sort descending by LZB (highest first)
  std::sort(topHotBuckets.begin(), topHotBuckets.end());

  // Keep only top 5
  if(topHotBuckets.size() > 5) {
    topHotBuckets.resize(5);
  }
}

// Calculate collision probability based on birthday paradox
// P(collision) = 1 - exp(-N^2 / (2 * M))
// where N = number of DPs, M = search space size
double Kangaroo::CalculateCollisionProbability(uint64_t numDPs) {

  if(numDPs == 0) return 0.0;

  // Use rangeWidth as search space M
  // For large numbers, use approximation to avoid overflow
  double N = (double)numDPs;
  double M = rangeWidth.ToDouble();

  if(M <= 0.0) return 0.0;

  // Calculate -N^2 / (2M)
  double exponent = -(N * N) / (2.0 * M);

  // For very small exponents (large negative), exp() approaches 0, so P approaches 1
  if(exponent < -100.0) return 1.0;

  // For very large exponents (close to 0), use approximation
  if(exponent > -0.001) {
    // P ≈ N^2 / (2M) for small probabilities
    return (N * N) / (2.0 * M);
  }

  // Normal case
  double prob = 1.0 - exp(exponent);
  return (prob < 0.0) ? 0.0 : ((prob > 1.0) ? 1.0 : prob);
}

// Calculate ETA to reach a target collision probability
// Returns time in seconds, or -1.0 if target already reached or invalid
double Kangaroo::CalculateETAForProbability(double targetProb, uint64_t currentDPs, double currentRate) {

  if(targetProb <= 0.0 || targetProb >= 1.0) return -1.0;
  if(currentRate <= 0.0) return -1.0;

  // Check if we've already reached the target
  double currentProb = CalculateCollisionProbability(currentDPs);
  if(currentProb >= targetProb) return 0.0;

  // Solve for N where P(N) = targetProb
  // P = 1 - exp(-N^2 / (2M))
  // targetProb = 1 - exp(-N^2 / (2M))
  // exp(-N^2 / (2M)) = 1 - targetProb
  // -N^2 / (2M) = ln(1 - targetProb)
  // N^2 = -2M * ln(1 - targetProb)
  // N = sqrt(-2M * ln(1 - targetProb))

  double M = rangeWidth.ToDouble();
  if(M <= 0.0) return -1.0;

  double logTerm = log(1.0 - targetProb);
  if(logTerm >= 0.0) return -1.0; // Invalid

  double targetN = sqrt(-2.0 * M * logTerm);

  if(targetN <= (double)currentDPs) return 0.0;

  // Calculate remaining DPs needed
  double remainingDPs = targetN - (double)currentDPs;

  // Calculate time based on current rate
  double timeSeconds = remainingDPs / currentRate;

  return timeSeconds;
}

// Wait for end of server and dispay stats
void Kangaroo::ProcessServer() {

  double t0;
  double t1;
  t0 = Timer::get_tick();
  startTime = t0;
  double lastSave = 0;

  // Acquire mutex ownership
#ifndef WIN64
  pthread_mutex_init(&ghMutex, NULL);
  setvbuf(stdout, NULL, _IONBF, 0);
#else
  ghMutex = CreateMutex(NULL,FALSE,NULL);
#endif

  while(!endOfSearch) {

    t0 = Timer::get_tick();

    LOCK(ghMutex);
    // Get back all dps
    localCache.clear();
    for(int i=0;i<(int)recvDP.size();i++)
      localCache.push_back(recvDP[i]);
    recvDP.clear();
    UNLOCK(ghMutex);

    // Add to hashTable
    for(int i = 0; i<(int)localCache.size() && !endOfSearch; i++) {
      DP_CACHE dp = localCache[i];
      for(int j = 0; j<(int)dp.nbDP && !endOfSearch; j++) {
        uint64_t h = dp.dp[j].h;
        if(!AddToTable(h,&dp.dp[j].x,&dp.dp[j].d)) {
          // Collision inside the same herd
          collisionInSameHerd++;
        }
      }
      free(dp.dp);
    }

    t1 = Timer::get_tick();

    double toSleep = SEND_PERIOD - (t1-t0);
    if(toSleep<0) toSleep = 0.0;
    Timer::SleepMillis((uint32_t)(toSleep*1000.0));

    t1 = Timer::get_tick();

    if(!endOfSearch) {
      // Calculate tame/wild ratio (1.000 = 50/50)
      double twRatio = (wildCount > 0) ? ((double)tameCount / (double)wildCount) : 0.0;

      // Calculate compact gap values - convert to billions
      // Full 128-bit value: (high * 2^64 + low) / 1e9
      double gap128 = (double)lastGap.i64[1] * 18446744073709551616.0 + (double)lastGap.i64[0];
      double lowestGap128 = (double)lowestGap.i64[1] * 18446744073709551616.0 + (double)lowestGap.i64[0];
      double currentGap = gap128 / 1000000000.0;
      double lowest = lowestGap128 / 1000000000.0;

      // First line - original progress bar
      printf("\r[Client %d][Kang 2^%.2f][DP Count 2^%.2f/2^%.2f][Dead %.0f][T/W:%.3f][Gap:%.1f][L.Gap:%.1f][%s][%s]  ",
        connectedClient,
        log2((double)totalRW),
        log2((double)hashTable.GetNbItem()),
        log2(expectedNbOp / pow(2.0,dpSize)),
        (double)collisionInSameHerd,
        twRatio,
        currentGap, lowest,
        GetTimeStr(t1 - startTime).c_str(),
        hashTable.GetSizeInfo().c_str()
        );

      // Second line - LZB analysis and collision probability
      uint64_t numDPs = hashTable.GetNbItem();
      double collisionProb = CalculateCollisionProbability(numDPs);

      // Calculate expected max LZB based on birthday paradox
      double expectedMaxLZB = (numDPs > 1) ? log2((double)numDPs * (double)numDPs / 2.0) : 0.0;

      // Calculate 50% and 90% ETA
      double dpRate = (t1 - startTime > 0) ? (double)numDPs / (t1 - startTime) : 0.0;
      double eta50 = CalculateETAForProbability(0.50, numDPs, dpRate);
      double eta90 = CalculateETAForProbability(0.90, numDPs, dpRate);

      // Format hot buckets display
      char hotBucketsStr[128] = "";
      LOCK(ghMutex);
      if(topHotBuckets.size() > 0) {
        char temp[32];
        strcpy(hotBucketsStr, "Hot:[");
        for(size_t i = 0; i < topHotBuckets.size() && i < 3; i++) {
          if(i > 0) strcat(hotBucketsStr, ",");
          sprintf(temp, "0x%X(%u)", topHotBuckets[i].bucketId, topHotBuckets[i].lzb);
          strcat(hotBucketsStr, temp);
        }
        strcat(hotBucketsStr, "]");
      }
      uint32_t localMaxLZB = maxLeadingZeroBits;
      UNLOCK(ghMutex);

      printf("\n\033[K[LZB:%u/%.0f][P:%.1f%%][50%%:%s][90%%:%s]%s  \033[F",
        localMaxLZB,
        expectedMaxLZB,
        collisionProb * 100.0,
        (eta50 >= 0) ? GetTimeStr(eta50).c_str() : "N/A",
        (eta90 >= 0) ? GetTimeStr(eta90).c_str() : "N/A",
        hotBucketsStr
        );
      fflush(stdout);
    }

    if(workFile.length() > 0 && !endOfSearch) {
      if((t1 - lastSave) > saveWorkPeriod) {
        SaveServerWork();
        lastSave = t1;
      }
    }

  }

}

// Wait for end of threads and display stats
void Kangaroo::Process(TH_PARAM *params,std::string unit) {

  double t0;
  double t1;

  uint64_t count;
  uint64_t lastCount = 0;
  uint64_t gpuCount = 0;
  uint64_t lastGPUCount = 0;
  double avgKeyRate = 0.0;
  double avgGpuKeyRate = 0.0;
  double lastSave = 0;

#ifndef WIN64
  setvbuf(stdout, NULL, _IONBF, 0);
#endif

  // Key rate smoothing filter
#define FILTER_SIZE 8
  double lastkeyRate[FILTER_SIZE];
  double lastGpukeyRate[FILTER_SIZE];
  uint32_t filterPos = 0;

  double keyRate = 0.0;
  double gpuKeyRate = 0.0;

  memset(lastkeyRate,0,sizeof(lastkeyRate));
  memset(lastGpukeyRate,0,sizeof(lastkeyRate));

  // Wait that all threads have started
  while(!hasStarted(params))
    Timer::SleepMillis(5);

  t0 = Timer::get_tick();
  startTime = t0;
  lastGPUCount = getGPUCount();
  lastCount = getCPUCount() + gpuCount;

  while(isAlive(params)) {

    int delay = 2000;
    while(isAlive(params) && delay > 0) {
      Timer::SleepMillis(50);
      delay -= 50;
    }

    gpuCount = getGPUCount();
    count = getCPUCount() + gpuCount;

    t1 = Timer::get_tick();
    keyRate = (double)(count - lastCount) / (t1 - t0);
    gpuKeyRate = (double)(gpuCount - lastGPUCount) / (t1 - t0);
    lastkeyRate[filterPos%FILTER_SIZE] = keyRate;
    lastGpukeyRate[filterPos%FILTER_SIZE] = gpuKeyRate;
    filterPos++;

    // KeyRate smoothing
    uint32_t nbSample;
    for(nbSample = 0; (nbSample < FILTER_SIZE) && (nbSample < filterPos); nbSample++) {
      avgKeyRate += lastkeyRate[nbSample];
      avgGpuKeyRate += lastGpukeyRate[nbSample];
    }
    avgKeyRate /= (double)(nbSample);
    avgGpuKeyRate /= (double)(nbSample);
    double expectedTime = expectedNbOp / avgKeyRate;

    // Display stats
    if(isAlive(params) && !endOfSearch) {
      // Calculate tame/wild ratio (1.000 = 50/50)
      double twRatio = (wildCount > 0) ? ((double)tameCount / (double)wildCount) : 0.0;

      // Calculate compact gap values - convert to billions
      // Full 128-bit value: (high * 2^64 + low) / 1e9
      double gap128 = (double)lastGap.i64[1] * 18446744073709551616.0 + (double)lastGap.i64[0];
      double lowestGap128 = (double)lowestGap.i64[1] * 18446744073709551616.0 + (double)lowestGap.i64[0];
      double currentGap = gap128 / 1000000000.0;
      double lowest = lowestGap128 / 1000000000.0;

      // First line - original progress bar
      if(clientMode) {
        printf("\r[%.2f %s][GPU %.2f %s][Count 2^%.2f][T/W:%.3f][Gap:%.1f][L.Gap:%.1f][%s][Server %6s]  ",
          avgKeyRate / 1000000.0,unit.c_str(),
          avgGpuKeyRate / 1000000.0,unit.c_str(),
          log2((double)count + offsetCount),
          twRatio,
          currentGap, lowest,
          GetTimeStr(t1 - startTime + offsetTime).c_str(),
          serverStatus.c_str()
          );
      } else {
        printf("\r[%.2f %s][GPU %.2f %s][Count 2^%.2f][Dead %.0f][T/W:%.3f][Gap:%.1f][L.Gap:%.1f][%s (Avg %s)][%s]  ",
          avgKeyRate / 1000000.0,unit.c_str(),
          avgGpuKeyRate / 1000000.0,unit.c_str(),
          log2((double)count + offsetCount),
          (double)collisionInSameHerd,
          twRatio,
          currentGap, lowest,
          GetTimeStr(t1 - startTime + offsetTime).c_str(),GetTimeStr(expectedTime).c_str(),
          hashTable.GetSizeInfo().c_str()
        );
      }

      // Second line - LZB analysis and collision probability
      uint64_t numDPs = hashTable.GetNbItem();
      double collisionProb = CalculateCollisionProbability(numDPs);

      // Calculate expected max LZB based on birthday paradox
      double expectedMaxLZB = (numDPs > 1) ? log2((double)numDPs * (double)numDPs / 2.0) : 0.0;

      // Calculate 50% and 90% ETA (use count rate, not DP rate for better estimation)
      double currentTime = t1 - startTime + offsetTime;
      double countRate = (currentTime > 0) ? ((double)count + offsetCount) / currentTime : 0.0;
      double dpRate = (currentTime > 0) ? (double)numDPs / currentTime : 0.0;
      double eta50 = CalculateETAForProbability(0.50, numDPs, dpRate);
      double eta90 = CalculateETAForProbability(0.90, numDPs, dpRate);

      // Format hot buckets display
      char hotBucketsStr[128] = "";
      LOCK(ghMutex);
      if(topHotBuckets.size() > 0) {
        char temp[32];
        strcpy(hotBucketsStr, "Hot:[");
        for(size_t i = 0; i < topHotBuckets.size() && i < 3; i++) {
          if(i > 0) strcat(hotBucketsStr, ",");
          sprintf(temp, "0x%X(%u)", topHotBuckets[i].bucketId, topHotBuckets[i].lzb);
          strcat(hotBucketsStr, temp);
        }
        strcat(hotBucketsStr, "]");
      }
      uint32_t localMaxLZB = maxLeadingZeroBits;
      UNLOCK(ghMutex);

      printf("\n\033[K[LZB:%u/%.0f][P:%.1f%%][50%%:%s][90%%:%s]%s  \033[F",
        localMaxLZB,
        expectedMaxLZB,
        collisionProb * 100.0,
        (eta50 >= 0) ? GetTimeStr(eta50).c_str() : "N/A",
        (eta90 >= 0) ? GetTimeStr(eta90).c_str() : "N/A",
        hotBucketsStr
        );
      fflush(stdout);

    }

    // Save request
    if(workFile.length() > 0 && !endOfSearch) {
      if((t1 - lastSave) > saveWorkPeriod) {
        SaveWork(count + offsetCount,t1 - startTime + offsetTime,params,nbCPUThread + nbGPUThread);
        lastSave = t1;
      }
    }

    // Abort
    if(!clientMode && maxStep>0.0) {
      double max = expectedNbOp * maxStep; 
      if( (double)count > max ) {
        ::printf("\nKey#%2d [XX]Pub:  0x%s \n",keyIdx,secp->GetPublicKeyHex(true,keysToSearch[keyIdx]).c_str());
        ::printf("       Aborted !\n");
        endOfSearch = true;
        Timer::SleepMillis(1000);
      }
    }

    lastCount = count;
    lastGPUCount = gpuCount;
    t0 = t1;

  }

  count = getCPUCount() + getGPUCount();
  t1 = Timer::get_tick();

  if( !endOfSearch ) {
    printf("\r[%.2f %s][GPU %.2f %s][Cnt 2^%.2f][%s]  ",
      avgKeyRate / 1000000.0,unit.c_str(),
      avgGpuKeyRate / 1000000.0,unit.c_str(),
      log2((double)count),
      GetTimeStr(t1 - startTime).c_str()
      );
  }

}

// ----------------------------------------------------------------------------

void Kangaroo::ScanGapsThread(TH_PARAM *p) {

  // Background thread for periodic gap tracking
  // This runs independently and doesn't affect the critical path

#ifndef WIN64
  setvbuf(stdout, NULL, _IONBF, 0);
#endif

  while(!endOfSearch) {

    // Sleep for 3 seconds between scans
    int delay = 3000;
    while(!endOfSearch && delay > 0) {
      Timer::SleepMillis(50);
      delay -= 50;
    }

    if(endOfSearch) break;

    // Scan hash table for gaps
    int128_t localMinGap;
    localMinGap.i64[0] = 0xFFFFFFFFFFFFFFFFULL;
    localMinGap.i64[1] = 0x3FFFFFFFFFFFFFFFULL;

    int128_t localLastGap = lastGap;
    bool gapFound = false;

    // LZB analysis tracking
    uint32_t localMaxLZB = 0;
    int128_t localClosestXor;
    localClosestXor.i64[0] = 0xFFFFFFFFFFFFFFFFULL;
    localClosestXor.i64[1] = 0xFFFFFFFFFFFFFFFFULL;
    uint32_t localMaxLZBBucket = 0;

    std::vector<int128_t> distances;
    std::vector<int128_t> positions;    // X-coordinates for LZB analysis
    std::vector<uint32_t> herdTypes;

    // Scan through all hash buckets
    for(uint32_t h = 0; h < HASH_SIZE && !endOfSearch; h++) {

      distances.clear();
      positions.clear();
      herdTypes.clear();

      LOCK(ghMutex);
      uint32_t nbItem = hashTable.E[h].nbItem;

      if(nbItem > 1) {
        distances.reserve(nbItem);
        positions.reserve(nbItem);
        herdTypes.reserve(nbItem);

        for(uint32_t i = 0; i < nbItem; i++) {
          ENTRY* entry = hashTable.E[h].items[i];
          uint32_t type = (entry->d.i64[1] & 0x4000000000000000ULL) != 0;
          int128_t dist = entry->d;
          dist.i64[1] &= 0x3FFFFFFFFFFFFFFFULL;

          distances.push_back(dist);
          positions.push_back(entry->x);  // Store x-coordinate
          herdTypes.push_back(type);
        }
      }
      UNLOCK(ghMutex);

      uint32_t nbStored = (uint32_t)distances.size();
      if(nbStored > 1 && !endOfSearch) {

        // Track best LZB for this bucket
        uint32_t bucketMaxLZB = 0;
        int128_t bucketBestXor;
        bucketBestXor.i64[0] = 0xFFFFFFFFFFFFFFFFULL;
        bucketBestXor.i64[1] = 0xFFFFFFFFFFFFFFFFULL;

        for(uint32_t i = 0; i < nbStored - 1 && !endOfSearch; i++) {
          for(uint32_t j = i + 1; j < nbStored && !endOfSearch; j++) {
            if(herdTypes[i] != herdTypes[j]) {
              // Calculate absolute difference in distances (existing gap calculation)
              int128_t gap;
              if(distances[i].i64[1] > distances[j].i64[1] ||
                 (distances[i].i64[1] == distances[j].i64[1] && distances[i].i64[0] > distances[j].i64[0])) {
                gap.i64[1] = distances[i].i64[1] - distances[j].i64[1];
                gap.i64[0] = distances[i].i64[0] - distances[j].i64[0];
                if(distances[i].i64[0] < distances[j].i64[0]) gap.i64[1]--;
              } else {
                gap.i64[1] = distances[j].i64[1] - distances[i].i64[1];
                gap.i64[0] = distances[j].i64[0] - distances[i].i64[0];
                if(distances[j].i64[0] < distances[i].i64[0]) gap.i64[1]--;
              }

              gapFound = true;
              localLastGap = gap;

              // Update local minimum gap
              if(gap.i64[1] < localMinGap.i64[1] ||
                 (gap.i64[1] == localMinGap.i64[1] && gap.i64[0] < localMinGap.i64[0])) {
                localMinGap.i64[0] = gap.i64[0];
                localMinGap.i64[1] = gap.i64[1];
              }

              // NEW: Calculate XOR of x-coordinates for LZB analysis
              int128_t xorDiff;
              xorDiff.i64[0] = positions[i].i64[0] ^ positions[j].i64[0];
              xorDiff.i64[1] = positions[i].i64[1] ^ positions[j].i64[1];

              // Count leading zero bits
              uint32_t lzb = CountLeadingZeroBits(&xorDiff);

              // Update bucket-level best
              if(lzb > bucketMaxLZB) {
                bucketMaxLZB = lzb;
                bucketBestXor = xorDiff;
              }

              // Update global best
              if(lzb > localMaxLZB) {
                localMaxLZB = lzb;
                localClosestXor = xorDiff;
                localMaxLZBBucket = h;
              }
            }
          }
        }

        // Update hot buckets if this bucket has meaningful LZB
        if(bucketMaxLZB > 0) {
          UpdateHotBuckets(h, bucketMaxLZB, &bucketBestXor);
        }
      }
    }

    // Update global minimum gap, lowest gap, and last seen gap
    if(gapFound) {
      LOCK(ghMutex);
      lastGap.i64[0] = localLastGap.i64[0];
      lastGap.i64[1] = localLastGap.i64[1];

      minGap.i64[0] = localMinGap.i64[0];
      minGap.i64[1] = localMinGap.i64[1];

      // Update lowestGap only if this is a new all-time minimum
      if(localMinGap.i64[1] < lowestGap.i64[1] ||
         (localMinGap.i64[1] == lowestGap.i64[1] && localMinGap.i64[0] < lowestGap.i64[0])) {
        lowestGap.i64[0] = localMinGap.i64[0];
        lowestGap.i64[1] = localMinGap.i64[1];
      }

      // Update LZB analysis results
      if(localMaxLZB > maxLeadingZeroBits) {
        maxLeadingZeroBits = localMaxLZB;
        closestXorDiff.i64[0] = localClosestXor.i64[0];
        closestXorDiff.i64[1] = localClosestXor.i64[1];
        maxLZBBucket = localMaxLZBBucket;
      }

      UNLOCK(ghMutex);
    }

  }

}

