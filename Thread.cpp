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

  uint32_t scanOffset = 0;
  uint32_t bucketsPerScan = HASH_SIZE >> 3; // scan 1/8th of the table per pass
  if(bucketsPerScan == 0) bucketsPerScan = 1;

  std::vector<int256_t> distances;
  std::vector<int256_t> rawDistances;
  std::vector<uint32_t> herdTypes;

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
        if(!AddToTable(h,&dp.dp[j].x,&dp.dp[j].d,dp.dp[j].kType)) {
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
      // Snapshot gap and key estimate information together to keep display values in sync
      int256_t displayLastGap;
      int256_t displayLowestGap;
      Int displayKeyEstimate;
      bool displayHasKey = false;

      LOCK(ghMutex);
      displayLastGap = lastGap;
      displayLowestGap = lowestGap;
      displayHasKey = hasKeyEstimate;
      if(displayHasKey) {
        displayKeyEstimate.Set(&lastKeyEstimate);
      }
      UNLOCK(ghMutex);

      // Calculate tame/wild ratio (1.000 = 50/50)
      double twRatio = (wildCount > 0) ? ((double)tameCount / (double)wildCount) : 0.0;

      // Calculate compact gap values - convert to billions using the full 256-bit value
      auto toBillions = [](const int256_t& gap) {
        long double v = (long double)gap.i64[3];
        v = v * 18446744073709551616.0L + (long double)gap.i64[2];
        v = v * 18446744073709551616.0L + (long double)gap.i64[1];
        v = v * 18446744073709551616.0L + (long double)gap.i64[0];
        return (double)(v / 1000000000.0L);
      };

      double currentGap = toBillions(displayLastGap);
      double lowest = toBillions(displayLowestGap);

      // Second line - DP insertion metrics
      uint64_t currentDPs = hashTable.GetNbItem();
      double currentTime = t1 - startTime;

      // Calculate DP rate (using smoothed calculation)
      double dpRate = 0.0;
      if(currentTime > 0.0) {
        // Calculate instantaneous rate based on change since last update
        if(lastDPTime > 0.0 && (currentTime - lastDPTime) > 0.0) {
          dpRate = (double)(currentDPs - lastDPCount) / (currentTime - lastDPTime);
        } else {
          // First time or no time elapsed, use average rate
          dpRate = (double)currentDPs / currentTime;
        }
      }

      // Update tracking variables for next iteration
      lastDPCount = currentDPs;
      lastDPTime = currentTime;

      std::string keyEstimateStr = displayHasKey ? displayKeyEstimate.GetBase10() : "n/a";

      // First line - original progress bar
      printf("\r[Client %d][Kang 2^%.2f][Dead %.0f][T/W:%.3f][L.Gap:%.1f][k_est:%s][%s]  ",
        connectedClient,
        log2((double)totalRW),
        (double)collisionInSameHerd,
        twRatio,
        lowest,
        keyEstimateStr.c_str(),
        GetTimeStr(t1 - startTime).c_str()
        );

      printf("\n\033[K[DP Count 2^%.2f/2^%.2f][Gap:%.1f][%s][DP Inserted: %llu][DP Rate: %.2f DP/s]%s  \033[F",
        log2((double)hashTable.GetNbItem()),
        log2(expectedNbOp / pow(2.0,dpSize)),
        currentGap,
        hashTable.GetSizeInfo().c_str(),
        (unsigned long long)currentDPs,
        dpRate,
        maxStepExtended ? "[extended Scan]" : ""
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

      auto toBillions = [](const int256_t& gap) {
        long double v = (long double)gap.i64[3];
        v = v * 18446744073709551616.0L + (long double)gap.i64[2];
        v = v * 18446744073709551616.0L + (long double)gap.i64[1];
        v = v * 18446744073709551616.0L + (long double)gap.i64[0];
        return (double)(v / 1000000000.0L);
      };

      double currentGap = toBillions(lastGap);
      double lowest = toBillions(lowestGap);

      // Second line - DP insertion metrics
      uint64_t currentDPs = hashTable.GetNbItem();
      double currentTime = t1 - startTime + offsetTime;

      // Calculate DP rate (using smoothed calculation)
      double dpRate = 0.0;
      if(currentTime > 0.0) {
        // Calculate instantaneous rate based on change since last update
        if(lastDPTime > 0.0 && (currentTime - lastDPTime) > 0.0) {
          dpRate = (double)(currentDPs - lastDPCount) / (currentTime - lastDPTime);
        } else {
          // First time or no time elapsed, use average rate
          dpRate = (double)currentDPs / currentTime;
        }
      }

      // Update tracking variables for next iteration
      lastDPCount = currentDPs;
      lastDPTime = currentTime;

      // Snapshot the key estimate alongside the latest gap so it matches the displayed L.Gap
      Int displayKeyEstimate;
      bool displayHasKey = false;

      LOCK(ghMutex);
      displayHasKey = hasKeyEstimate;
      if(displayHasKey) {
        displayKeyEstimate.Set(&lastKeyEstimate);
      }
      UNLOCK(ghMutex);

      std::string keyEstimateStr = displayHasKey ? displayKeyEstimate.GetBase10() : "n/a";

      // First line - original progress bar
      if(clientMode) {
        printf("\r[%.2f %s][GPU %.2f %s][T/W:%.3f][L.Gap:%.1f][k_est:%s][%s][Server %6s]  ",
          avgKeyRate / 1000000.0,unit.c_str(),
          avgGpuKeyRate / 1000000.0,unit.c_str(),
          twRatio,
          lowest,
          keyEstimateStr.c_str(),
          GetTimeStr(t1 - startTime + offsetTime).c_str(),
          serverStatus.c_str()
          );
      } else {
        printf("\r[%.2f %s][GPU %.2f %s][Dead %.0f][T/W:%.3f][L.Gap:%.1f][k_est:%s][%s (Avg %s)]  ",
          avgKeyRate / 1000000.0,unit.c_str(),
          avgGpuKeyRate / 1000000.0,unit.c_str(),
          (double)collisionInSameHerd,
          twRatio,
          lowest,
          keyEstimateStr.c_str(),
          GetTimeStr(t1 - startTime + offsetTime).c_str(),GetTimeStr(expectedTime).c_str()
        );
      }

      printf("\n\033[K[Count 2^%.2f][Gap:%.1f][%s][DP Inserted: %llu][DP Rate: %.2f DP/s]%s  \033[F",
        log2((double)count + offsetCount),
        currentGap,
        hashTable.GetSizeInfo().c_str(),
        (unsigned long long)currentDPs,
        dpRate,
        maxStepExtended ? "[extended Scan]" : ""
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
        // Check if we can extend based on l.gap
        if(!maxStepExtended) {
          // Calculate lowestGap in decimal form (divide by 1e9) using all 256 bits
          long double lowestGapValue = (long double)lowestGap.i64[3];
          lowestGapValue = lowestGapValue * 18446744073709551616.0L + (long double)lowestGap.i64[2];
          lowestGapValue = lowestGapValue * 18446744073709551616.0L + (long double)lowestGap.i64[1];
          lowestGapValue = lowestGapValue * 18446744073709551616.0L + (long double)lowestGap.i64[0];
          double lowestGapDecimal = (double)(lowestGapValue / 1000000000.0L);

          if(lowestGapDecimal < 1000.0 && lowestGapDecimal > 0.0) {
            // Extend scan by +2
            maxStep += 2.0;
            maxStepExtended = true;
            // Continue scanning - don't abort yet
          } else {
            // Gap too large or not available, abort
            ::printf("\n\nKey#%2d [XX]Pub:  0x%s \n",keyIdx,secp->GetPublicKeyHex(true,keysToSearch[keyIdx]).c_str());
            ::printf("       Aborted !\n");
            endOfSearch = true;
            Timer::SleepMillis(1000);
          }
        } else {
          // Already extended once, abort now
          ::printf("\n\nKey#%2d [XX]Pub:  0x%s \n",keyIdx,secp->GetPublicKeyHex(true,keysToSearch[keyIdx]).c_str());
          ::printf("       Aborted !\n");
          endOfSearch = true;
          Timer::SleepMillis(1000);
        }
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

  // Reuse buffers across scans to avoid repeated allocations
  static uint32_t scanOffset = 0;
  static uint32_t bucketsPerScan = (HASH_SIZE >> 3) ? (HASH_SIZE >> 3) : 1; // scan 1/8th of the table per pass
  static std::vector<int256_t> distances;
  static std::vector<int256_t> rawDistances;
  static std::vector<uint32_t> herdTypes;

  while(!endOfSearch) {

    // Sleep for 3 seconds between scans
    int delay = 3000;
    while(!endOfSearch && delay > 0) {
      Timer::SleepMillis(50);
      delay -= 50;
    }

    if(endOfSearch) break;

    // Scan hash table for gaps
    int256_t localMinGap;
    localMinGap.i64[0] = 0xFFFFFFFFFFFFFFFFULL;
    localMinGap.i64[1] = 0xFFFFFFFFFFFFFFFFULL;
    localMinGap.i64[2] = 0xFFFFFFFFFFFFFFFFULL;
    localMinGap.i64[3] = 0xFFFFFFFFFFFFFFFFULL;

    int256_t localLastGap = lastGap;
    bool gapFound = false;
    Int localKeyEstimate;
    localKeyEstimate.SetInt32(0);

    // Track the key estimate associated with the smallest gap found in this scan
    Int minGapKeyEstimate;
    minGapKeyEstimate.SetInt32(0);
    bool hasMinGapKey = false;

    // Scan only a subset of the hash buckets per pass to reduce contention
    uint32_t startBucket = scanOffset;
    uint32_t endBucket = startBucket + bucketsPerScan;
    for(uint32_t h = startBucket; h < endBucket && !endOfSearch; h++) {
      uint32_t bucket = (h & HASH_MASK);

      distances.clear();
      rawDistances.clear();
      herdTypes.clear();

      LOCK(ghMutex);
      uint32_t nbItem = hashTable.E[bucket].nbItem;

      if(nbItem > 1) {
        distances.reserve(nbItem);
        rawDistances.reserve(nbItem);
        herdTypes.reserve(nbItem);

        for(uint32_t i = 0; i < nbItem; i++) {
          ENTRY* entry = hashTable.E[bucket].items[i];
          int256_t dist = entry->d;
          uint32_t type = entry->kType;
          rawDistances.push_back(dist);

          distances.push_back(dist);
          herdTypes.push_back(type);
        }
      }
      UNLOCK(ghMutex);

      uint32_t nbStored = (uint32_t)distances.size();
      if(nbStored > 1 && !endOfSearch) {

        for(uint32_t i = 0; i < nbStored - 1 && !endOfSearch; i++) {
          for(uint32_t j = i + 1; j < nbStored && !endOfSearch; j++) {
            if(herdTypes[i] != herdTypes[j]) {
              uint32_t tameIdx = (herdTypes[i] == TAME) ? i : j;
              uint32_t wildIdx = (herdTypes[i] == WILD) ? i : j;

              // Calculate absolute difference in distances using Int math so
              // the gap follows the same modulo semantics as collision checks
              Int tameDistance;
              Int wildDistance;
              HashTable::CalcDist(&rawDistances[tameIdx], &tameDistance);
              HashTable::CalcDist(&rawDistances[wildIdx], &wildDistance);

              // Calculate gap using modular arithmetic to handle wraparound correctly
              // The gap is the minimum distance between two points on the modular circle
              Int gapInt;
              gapInt.Set(&tameDistance);
              gapInt.ModSubK1order(&wildDistance);  // Compute (tame - wild) mod order

              // Check if the opposite direction is shorter
              // If gap > order/2, use (order - gap) instead
              Int gapNeg;
              gapNeg.Set(&gapInt);
              gapNeg.ModNegK1order();  // gapNeg = order - gap (modular negation)

              // Use the smaller of the two distances
              if(gapNeg.IsLower(&gapInt)) {
                gapInt.Set(&gapNeg);
              }

              int256_t gap;
              gap.i64[0] = gapInt.bits64[0];
              gap.i64[1] = gapInt.bits64[1];
              gap.i64[2] = gapInt.bits64[2];
              gap.i64[3] = gapInt.bits64[3];

              gapFound = true;
              localLastGap = gap;

              // Compute key estimate using the same approach as CheckKey()
              // Try all four sign combinations (type 0-3) when adding the two distances
              bool hasKeyCandidate = false;
              for(uint8_t type = 0; type < 4; type++) {

                Int tdAdj(&tameDistance);
                Int wdAdj(&wildDistance);

                if(type & 0x1)
                  tdAdj.ModNegK1order();
                if(type & 0x2)
                  wdAdj.ModNegK1order();

                Int keyEst(&tdAdj);
                keyEst.ModAddK1order(&wdAdj);
#ifdef USE_SYMMETRY
                keyEst.ModAddK1order(&rangeWidthDiv2);
#endif
                keyEst.ModAddK1order(&rangeStart);

                bool in_range = !keyEst.IsLower(&rangeStart) && !keyEst.IsGreater(&rangeEnd);

                if(in_range || !hasKeyCandidate) {
                  localKeyEstimate.Set(&keyEst);
                  hasKeyCandidate = true;
                }

                if(in_range) {
                  break;
                }
              }

              // Update local minimum gap and remember its key estimate
              bool is_smaller = false;
              if(gap.i64[3] != localMinGap.i64[3]) {
                is_smaller = gap.i64[3] < localMinGap.i64[3];
              } else if(gap.i64[2] != localMinGap.i64[2]) {
                is_smaller = gap.i64[2] < localMinGap.i64[2];
              } else if(gap.i64[1] != localMinGap.i64[1]) {
                is_smaller = gap.i64[1] < localMinGap.i64[1];
              } else {
                is_smaller = gap.i64[0] < localMinGap.i64[0];
              }

              if(is_smaller) {
                localMinGap.i64[0] = gap.i64[0];
                localMinGap.i64[1] = gap.i64[1];
                localMinGap.i64[2] = gap.i64[2];
                localMinGap.i64[3] = gap.i64[3];

                if(hasKeyCandidate) {
                  minGapKeyEstimate.Set(&localKeyEstimate);
                  hasMinGapKey = true;
                }
              }
            }
          }
        }
      }
    }

    scanOffset = (scanOffset + bucketsPerScan) & HASH_MASK;

    // Update global minimum gap, lowest gap, and last seen gap
    if(gapFound) {
      LOCK(ghMutex);
      lastGap.i64[0] = localLastGap.i64[0];
      lastGap.i64[1] = localLastGap.i64[1];
      lastGap.i64[2] = localLastGap.i64[2];
      lastGap.i64[3] = localLastGap.i64[3];

      minGap.i64[0] = localMinGap.i64[0];
      minGap.i64[1] = localMinGap.i64[1];
      minGap.i64[2] = localMinGap.i64[2];
      minGap.i64[3] = localMinGap.i64[3];

      // Update lowestGap only if this is a new all-time minimum
      bool is_new_lowest = false;
      if(localMinGap.i64[3] != lowestGap.i64[3]) {
        is_new_lowest = localMinGap.i64[3] < lowestGap.i64[3];
      } else if(localMinGap.i64[2] != lowestGap.i64[2]) {
        is_new_lowest = localMinGap.i64[2] < lowestGap.i64[2];
      } else if(localMinGap.i64[1] != lowestGap.i64[1]) {
        is_new_lowest = localMinGap.i64[1] < lowestGap.i64[1];
      } else {
        is_new_lowest = localMinGap.i64[0] < lowestGap.i64[0];
      }

      if(is_new_lowest) {
        lowestGap.i64[0] = localMinGap.i64[0];
        lowestGap.i64[1] = localMinGap.i64[1];
        lowestGap.i64[2] = localMinGap.i64[2];
        lowestGap.i64[3] = localMinGap.i64[3];

        if(hasMinGapKey) {
          lastKeyEstimate.Set(&minGapKeyEstimate);
          hasKeyEstimate = true;
        }
      }

      UNLOCK(ghMutex);
    }

  }

}

