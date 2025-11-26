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

// ============================================================================
// GRADUATED DP STRATEGY IMPLEMENTATION
// ============================================================================

// Initialize graduated DP strategy
void Kangaroo::InitGraduatedDP(double actualKeyRate) {

  if(!gradConfig.enabled) {
    currentPhase = PHASE_DISABLED;
    return;
  }

  // Allocate per-bucket statistics
  if(bucketStats == nullptr) {
    bucketStats = new BucketStats[HASH_SIZE];
  }

  // Reset all bucket stats
  for(uint32_t i = 0; i < HASH_SIZE; i++) {
    bucketStats[i] = BucketStats();
  }

  // Reserve space for hotspots
  topHotspots.reserve(gradConfig.topHotspotsCount);

  // Calculate phase end times - use manual duration if specified, otherwise use actual performance
  double totalDuration;  // in seconds

  if(gradConfig.manualDuration > 0.0) {
    // User specified manual duration in hours
    totalDuration = gradConfig.manualDuration * 3600.0;  // Convert to seconds
    ::printf("\nUsing MANUAL duration: %.1f hours\n", gradConfig.manualDuration);
  } else {
    // Use actual measured key rate to calculate expected time
    if(actualKeyRate > 0.0) {
      totalDuration = expectedNbOp / actualKeyRate;
      ::printf("\nUsing ESTIMATED duration based on measured performance: %.2f hours\n", totalDuration / 3600.0);
    } else {
      // Fallback if no key rate available yet (shouldn't happen)
      totalDuration = expectedNbOp / 1e6;
      ::printf("\nUsing ESTIMATED duration (no performance data yet): %.2f hours\n", totalDuration / 3600.0);
    }
  }

  phase1EndTime = startTime + (totalDuration * gradConfig.phase1Duration);
  phase2EndTime = phase1EndTime + (totalDuration * gradConfig.phase2Duration);
  phase3EndTime = phase2EndTime + (totalDuration * gradConfig.phase3Duration);

  phaseStartTime = startTime;
  currentPhase = PHASE_WIDE_NET;

  ::printf("\n");
  ::printf("╔════════════════════════════════════════════════════════════════════════╗\n");
  ::printf("║           GRADUATED DP STRATEGY - MULTI-PHASE SEARCH ENABLED          ║\n");
  ::printf("╠════════════════════════════════════════════════════════════════════════╣\n");
  ::printf("║  Total Duration: %-52s     ║\n",
    gradConfig.manualDuration > 0.0 ?
      (std::string(GetTimeStr(totalDuration)) + " (manual)").c_str() :
      (std::string(GetTimeStr(totalDuration)) + " (estimated)").c_str());
  ::printf("╠════════════════════════════════════════════════════════════════════════╣\n");

  // Calculate actual DP bit counts for each phase with clamping
  int32_t phase1DP = (int32_t)dpSize + gradConfig.phase1DPBits;
  int32_t phase2DP = (int32_t)dpSize + gradConfig.phase2DPBits;
  int32_t phase3DP = (int32_t)dpSize + gradConfig.phase3DPBits;

  // Clamp to reasonable range [8, 64]
  if(phase1DP < 8) phase1DP = 8;
  if(phase1DP > 64) phase1DP = 64;
  if(phase2DP < 8) phase2DP = 8;
  if(phase2DP > 64) phase2DP = 64;
  if(phase3DP < 8) phase3DP = 8;
  if(phase3DP > 64) phase3DP = 64;

  ::printf("║  Phase 1: WIDE NET      - %.0f%% of time - DP: %d bits (FAST)          ║\n",
    gradConfig.phase1Duration * 100, phase1DP);
  ::printf("║           Duration: %-49s ║\n",
    GetTimeStr(totalDuration * gradConfig.phase1Duration).c_str());
  ::printf("║  Phase 2: FOCUSED       - %.0f%% of time - DP: %d bits (ADAPTIVE)      ║\n",
    gradConfig.phase2Duration * 100, phase2DP);
  ::printf("║           Duration: %-49s ║\n",
    GetTimeStr(totalDuration * gradConfig.phase2Duration).c_str());
  ::printf("║  Phase 3: PRECISION     - %.0f%% of time - DP: %d bits (TARGETED)      ║\n",
    gradConfig.phase3Duration * 100, phase3DP);
  ::printf("║           Duration: %-49s ║\n",
    GetTimeStr(totalDuration * gradConfig.phase3Duration).c_str());
  ::printf("╠════════════════════════════════════════════════════════════════════════╣\n");
  ::printf("║  Hotspot Bias Phase 2: %.0f%%  |  Hotspot Bias Phase 3: %.0f%%            ║\n",
    gradConfig.hotspotBiasPhase2 * 100, gradConfig.hotspotBiasPhase3 * 100);
  ::printf("║  Tracking Top %u Hotspots with %.1f%% decay per second               ║\n",
    gradConfig.topHotspotsCount, gradConfig.scoreDecayRate * 100);
  ::printf("╚════════════════════════════════════════════════════════════════════════╝\n");
  ::printf("\n");

  // Apply Phase 1 DP mask immediately
  uint32_t phase1DPSize = GetCurrentDPSize();
  SetDP(phase1DPSize);

  ResetPhaseStatistics();
}

// Update current phase based on time elapsed
void Kangaroo::UpdatePhase(double currentTime) {

  if(!gradConfig.enabled) return;

  SearchPhase oldPhase = currentPhase;

  if(currentTime >= phase3EndTime && !phase3Completed) {
    // Should not reach here, search should complete before phase 3 ends
    currentPhase = PHASE_PRECISION;
  } else if(currentTime >= phase2EndTime && !phase2Completed) {
    currentPhase = PHASE_PRECISION;
    if(oldPhase != PHASE_PRECISION) {
      phase2Completed = true;

      // Apply new DP mask for Phase 3
      uint32_t newDPSize = GetCurrentDPSize();
      SetDP(newDPSize);

      ::printf("\n");
      ::printf("╔════════════════════════════════════════════════════════════════════════╗\n");
      ::printf("║  ⚡ PHASE TRANSITION: FOCUSED → PRECISION                              ║\n");
      ::printf("║  Entering Phase 3: Precision Strike on hottest %u buckets            ║\n",
        (topHotspots.size() < 20) ? (uint32_t)topHotspots.size() : 20);
      ::printf("║  DP Mask: %u bits (ultra-dense coverage in hot zones)                 ║\n",
        newDPSize);
      ::printf("╚════════════════════════════════════════════════════════════════════════╝\n");
      ::printf("\n");
      ResetPhaseStatistics();
    }
  } else if(currentTime >= phase1EndTime && !phase1Completed) {
    currentPhase = PHASE_FOCUSED;
    if(oldPhase != PHASE_FOCUSED) {
      phase1Completed = true;

      // Apply new DP mask for Phase 2
      uint32_t newDPSize = GetCurrentDPSize();
      SetDP(newDPSize);

      ::printf("\n");
      ::printf("╔════════════════════════════════════════════════════════════════════════╗\n");
      ::printf("║  ⚡ PHASE TRANSITION: WIDE NET → FOCUSED                               ║\n");
      ::printf("║  Phase 1 Complete! Found %llu DPs in %.0f%% of expected time          ║\n",
        (unsigned long long)phase1DPCount, GetPhaseProgress() * 100);
      ::printf("║  Identified %u hotspots for focused search                            ║\n",
        (uint32_t)topHotspots.size());
      ::printf("║  DP Mask: %u bits with %.0f%% spawning in hotspots                     ║\n",
        newDPSize, gradConfig.hotspotBiasPhase2 * 100);
      ::printf("╚════════════════════════════════════════════════════════════════════════╝\n");
      ::printf("\n");
      ResetPhaseStatistics();
    }
  }

  // Update DP size if phase changed
  if(oldPhase != currentPhase) {
    uint32_t newDPSize = GetCurrentDPSize();
    SetDP(newDPSize);
  }
}

// Get current DP size based on phase
uint32_t Kangaroo::GetCurrentDPSize() {

  if(!gradConfig.enabled) {
    return initDPSize;  // Use original DP size
  }

  int32_t adjustment = 0;

  switch(currentPhase) {
    case PHASE_WIDE_NET:
      adjustment = gradConfig.phase1DPBits;
      break;
    case PHASE_FOCUSED:
      adjustment = gradConfig.phase2DPBits;
      break;
    case PHASE_PRECISION:
      adjustment = gradConfig.phase3DPBits;
      break;
    case PHASE_DISABLED:
    default:
      return initDPSize;
  }

  int32_t newSize = (int32_t)initDPSize + adjustment;

  // Clamp to reasonable range [8, 64]
  if(newSize < 8) newSize = 8;
  if(newSize > 64) newSize = 64;

  return (uint32_t)newSize;
}

// Update bucket statistics (called from ScanGapsThread)
void Kangaroo::UpdateBucketStatistics(double currentTime) {

  if(!gradConfig.enabled || bucketStats == nullptr) return;

  // This is called with ghMutex already locked from ScanGapsThread

  // Update stats for each bucket based on current hash table contents
  for(uint32_t h = 0; h < HASH_SIZE; h++) {

    uint32_t nbItem = hashTable.E[h].nbItem;
    if(nbItem == 0) continue;

    uint32_t tameInBucket = 0;
    uint32_t wildInBucket = 0;

    // Count TAME and WILD DPs in this bucket
    for(uint32_t i = 0; i < nbItem; i++) {
      ENTRY* entry = hashTable.E[h].items[i];
      uint32_t type = (entry->d.i64[1] & 0x4000000000000000ULL) != 0;
      if(type == 0) {
        tameInBucket++;
      } else {
        wildInBucket++;
      }
    }

    // Update counts
    bucketStats[h].tameCount = tameInBucket;
    bucketStats[h].wildCount = wildInBucket;

    // Update last arrival time if this bucket changed
    if(tameInBucket + wildInBucket > 0) {
      bucketStats[h].lastArrivalTime = currentTime;
    }

    // Calculate time-based decay factor (recent activity = higher weight)
    double timeSinceLastArrival = currentTime - bucketStats[h].lastArrivalTime;
    bucketStats[h].decayFactor = exp(-gradConfig.scoreDecayRate * timeSinceLastArrival);
  }
}

// Calculate hotspot scores for all buckets
void Kangaroo::CalculateHotspotScores() {

  if(!gradConfig.enabled || bucketStats == nullptr) return;

  uint64_t totalDPs = hashTable.GetNbItem();
  if(totalDPs == 0) return;

  double avgDPsPerBucket = (double)totalDPs / (double)HASH_SIZE;

  for(uint32_t h = 0; h < HASH_SIZE; h++) {

    uint32_t total = bucketStats[h].tameCount + bucketStats[h].wildCount;
    if(total == 0) {
      bucketStats[h].densityScore = 0.0;
      bucketStats[h].balanceScore = 0.0;
      bucketStats[h].collisionScore = 0.0;
      continue;
    }

    // Density score: how many DPs relative to average
    bucketStats[h].densityScore = (double)total / avgDPsPerBucket;

    // Balance score: how close is TAME/WILD ratio to 1.0
    if(bucketStats[h].tameCount > 0 && bucketStats[h].wildCount > 0) {
      double ratio = (double)bucketStats[h].tameCount / (double)bucketStats[h].wildCount;
      if(ratio > 1.0) ratio = 1.0 / ratio;  // Normalize to [0,1]
      bucketStats[h].balanceScore = ratio;  // 1.0 = perfect balance
    } else {
      bucketStats[h].balanceScore = 0.0;  // Only one herd present
    }

    // Collision score: combines density, balance, and recency
    // High score = high density + balanced TAME/WILD + recent activity
    bucketStats[h].collisionScore =
      bucketStats[h].densityScore *
      bucketStats[h].balanceScore *
      bucketStats[h].decayFactor;
  }
}

// Update top hotspots list
void Kangaroo::UpdateTopHotspots() {

  if(!gradConfig.enabled || bucketStats == nullptr) return;

  topHotspots.clear();

  // Collect all buckets with score above minimum threshold
  for(uint32_t h = 0; h < HASH_SIZE; h++) {
    if(bucketStats[h].collisionScore >= gradConfig.minHotspotScore) {
      Hotspot hs;
      hs.bucketId = h;
      hs.score = bucketStats[h].collisionScore;
      hs.tameCount = bucketStats[h].tameCount;
      hs.wildCount = bucketStats[h].wildCount;
      topHotspots.push_back(hs);
    }
  }

  // Sort by score (descending)
  std::sort(topHotspots.begin(), topHotspots.end());

  // Keep only top N
  if(topHotspots.size() > gradConfig.topHotspotsCount) {
    topHotspots.resize(gradConfig.topHotspotsCount);
  }
}

// Select spawn bucket based on current phase and hotspots
uint32_t Kangaroo::SelectSpawnBucket(bool isTame) {

  if(!gradConfig.enabled || currentPhase == PHASE_WIDE_NET || topHotspots.empty()) {
    // Phase 1 or no hotspots: uniform random spawning
    return (uint32_t)(rndl() % HASH_SIZE);
  }

  double bias = 0.0;
  if(currentPhase == PHASE_FOCUSED) {
    bias = gradConfig.hotspotBiasPhase2;
  } else if(currentPhase == PHASE_PRECISION) {
    bias = gradConfig.hotspotBiasPhase3;
  }

  // Decide whether to spawn in hotspot or uniformly
  double r = (double)rndl() / (double)0xFFFFFFFFFFFFFFFFULL;

  if(r < bias && !topHotspots.empty()) {
    // Spawn in a hotspot (weighted by score)

    // Calculate total score
    double totalScore = 0.0;
    for(size_t i = 0; i < topHotspots.size(); i++) {
      totalScore += topHotspots[i].score;
    }

    if(totalScore <= 0.0) {
      return (uint32_t)(rndl() % HASH_SIZE);  // Fallback to uniform
    }

    // Select hotspot weighted by score
    double target = ((double)rndl() / (double)0xFFFFFFFFFFFFFFFFULL) * totalScore;
    double cumulative = 0.0;

    for(size_t i = 0; i < topHotspots.size(); i++) {
      cumulative += topHotspots[i].score;
      if(cumulative >= target) {
        return topHotspots[i].bucketId;
      }
    }

    // Fallback (should not reach here)
    return topHotspots[0].bucketId;

  } else {
    // Uniform random spawning
    return (uint32_t)(rndl() % HASH_SIZE);
  }
}

// Respawn kangaroos to hotspot regions (ADAPTIVE RESEEDING - THE BRAIN)
void Kangaroo::RespawnKangaroosToHotspots(double respawnPercentage, TH_PARAM *threads, int nbThread) {

  if(!gradConfig.enabled || topHotspots.empty()) return;

  ::printf("║  Reseeding %.0f%% of kangaroos to %u hotspots...                          ║\n",
    respawnPercentage * 100, (uint32_t)topHotspots.size());

  uint64_t totalRespawned = 0;

  // Respawn kangaroos in each CPU thread
  for(int t = 0; t < nbThread; t++) {
    if(threads[t].isRunning) {
      uint64_t nbToRespawn = (uint64_t)(threads[t].nbKangaroo * respawnPercentage);

      for(uint64_t i = 0; i < nbToRespawn; i++) {
        // Select a random kangaroo to respawn
        uint64_t kangIdx = rndl() % threads[t].nbKangaroo;

        // Determine if this is TAME or WILD
        bool isTame = ((kangIdx + threads[t].threadId) % 2) == 0;

        // Select a hotspot bucket to spawn near
        uint32_t targetBucket = SelectSpawnBucket(isTame);

        // Generate new starting position near the target bucket
        // The bucket ID gives us a hint about the region of the search space
        Int newDistance;
        newDistance.Rand(rangePower);  // Random distance in search range

        // Bias the distance towards the hot zone
        // We use the bucket ID to add some spatial correlation
        Int bucketBias;
        bucketBias.SetInt32((int32_t)targetBucket);
        bucketBias.ShiftL(rangePower - 20);  // Scale appropriately
        newDistance.Add(&bucketBias);
        newDistance.Mod(&rangeWidth);

        // Set the new position
        threads[t].distance[kangIdx].Set(&newDistance);

        // Calculate the new point position
        Point newPoint;
        if(isTame) {
          // TAME: start + distance
          newPoint = secp->ComputePublicKey(&newDistance);
          newPoint = secp->AddDirect(keyToSearch, newPoint);
        } else {
          // WILD: distance only
          newPoint = secp->ComputePublicKey(&newDistance);
        }

        threads[t].px[kangIdx].Set(&newPoint.x);
        threads[t].py[kangIdx].Set(&newPoint.y);

        totalRespawned++;
      }
    }
  }

  ::printf("║  ✓ Respawned %llu kangaroos to hotspot regions                             ║\n",
    (unsigned long long)totalRespawned);
}

// Get phase information for display
void Kangaroo::GetPhaseInfo(char *buffer, size_t bufSize) {

  if(!gradConfig.enabled) {
    snprintf(buffer, bufSize, "STANDARD");
    return;
  }

  const char *phaseName = GetPhaseName(currentPhase);
  double progress = GetPhaseProgress() * 100.0;

  snprintf(buffer, bufSize, "%s[%.0f%%]", phaseName, progress);
}

// Reset phase statistics
void Kangaroo::ResetPhaseStatistics() {
  totalDPsAtPhaseStart = hashTable.GetNbItem();
}

// Get progress through current phase (0.0 to 1.0)
double Kangaroo::GetPhaseProgress() {

  if(!gradConfig.enabled) return 0.0;

  double currentTime = Timer::get_tick();
  double phaseStart, phaseEnd;

  switch(currentPhase) {
    case PHASE_WIDE_NET:
      phaseStart = startTime;
      phaseEnd = phase1EndTime;
      break;
    case PHASE_FOCUSED:
      phaseStart = phase1EndTime;
      phaseEnd = phase2EndTime;
      break;
    case PHASE_PRECISION:
      phaseStart = phase2EndTime;
      phaseEnd = phase3EndTime;
      break;
    default:
      return 0.0;
  }

  if(phaseEnd <= phaseStart) return 0.0;

  double progress = (currentTime - phaseStart) / (phaseEnd - phaseStart);
  if(progress < 0.0) progress = 0.0;
  if(progress > 1.0) progress = 1.0;

  return progress;
}

// Get phase name string
const char* Kangaroo::GetPhaseName(SearchPhase phase) {
  switch(phase) {
    case PHASE_WIDE_NET:  return "WIDE NET";
    case PHASE_FOCUSED:   return "FOCUSED";
    case PHASE_PRECISION: return "PRECISION";
    default:              return "STANDARD";
  }
}

// ============================================================================
// END GRADUATED DP STRATEGY
// ============================================================================

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

    // Initialize Graduated DP Strategy for server mode
    // Server mode should use manual duration (-gdp-time flag)
    if(gradConfig.enabled && !gradDPInitialized && !endOfSearch) {
      if((t1 - startTime) >= 10.0) {
        if(gradConfig.manualDuration > 0.0) {
          InitGraduatedDP(0.0);  // Server mode uses manual duration only
          gradDPInitialized = true;
        } else {
          // No manual duration specified in server mode, disable GDP
          ::printf("\nWARNING: Graduated DP requires -gdp-time flag in server mode. Disabling.\n");
          gradConfig.enabled = false;
        }
      }
    }

    // Update Graduated DP Strategy phases and bucket statistics
    if(gradConfig.enabled && gradDPInitialized && !endOfSearch) {
      UpdatePhase(t1);  // Check for phase transitions

      // Update bucket stats every 5 seconds
      if((t1 - lastHotspotUpdate) > 5.0) {
        UpdateBucketStatistics(t1);
        CalculateHotspotScores();
        UpdateTopHotspots();
        lastHotspotUpdate = t1;
      }
    }

    if(!endOfSearch) {
      // Calculate tame/wild ratio (1.000 = 50/50)
      double twRatio = (wildCount > 0) ? ((double)tameCount / (double)wildCount) : 0.0;

      // Calculate compact gap values - convert to billions
      // Full 128-bit value: (high * 2^64 + low) / 1e9
      double gap128 = (double)lastGap.i64[1] * 18446744073709551616.0 + (double)lastGap.i64[0];
      double lowestGap128 = (double)lowestGap.i64[1] * 18446744073709551616.0 + (double)lowestGap.i64[0];
      double currentGap = gap128 / 1000000000.0;
      // Check if lowestGap is still at initial max value (not yet set)
      bool lowestGapValid = !(lowestGap.i64[0] == 0xFFFFFFFFFFFFFFFFULL &&
                              lowestGap.i64[1] == 0x3FFFFFFFFFFFFFFFULL);
      double lowest = lowestGapValid ? (lowestGap128 / 1000000000.0) : 0.0;

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

      printf("\n\033[K[LZB:%u/%.0f][P:%.1f%%][50%%:%s][90%%:%s]%s  ",
        localMaxLZB,
        expectedMaxLZB,
        collisionProb * 100.0,
        (eta50 >= 0) ? GetTimeStr(eta50).c_str() : "N/A",
        (eta90 >= 0) ? GetTimeStr(eta90).c_str() : "N/A",
        hotBucketsStr
        );

      // Third line - Graduated DP phase info
      if(gradConfig.enabled) {
        char phaseInfo[128];
        GetPhaseInfo(phaseInfo, sizeof(phaseInfo));
        uint32_t topHotspotCount = (uint32_t)topHotspots.size();
        uint32_t topBucketId = topHotspots.empty() ? 0 : topHotspots[0].bucketId;
        double topScore = topHotspots.empty() ? 0.0 : topHotspots[0].score;

        printf("\n\033[K[PHASE:%s][Hotspots:%u][Top:0x%X(%.1f)]  \033[F\033[F",
          phaseInfo,
          topHotspotCount,
          topBucketId,
          topScore
          );
      } else {
        // Just move cursor back up if graduated DP is disabled
        printf("\033[F");
      }

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

    // Initialize Graduated DP Strategy once we have stable key rate measurements
    // Wait for at least 10 seconds and 3 samples for accurate estimation
    if(gradConfig.enabled && !gradDPInitialized && !endOfSearch) {
      if((t1 - startTime) >= 10.0 && nbSample >= 3 && avgKeyRate > 0.0) {
        InitGraduatedDP(avgKeyRate);
        gradDPInitialized = true;
      }
    }

    // Update Graduated DP Strategy phases and bucket statistics
    if(gradConfig.enabled && gradDPInitialized && !endOfSearch) {
      UpdatePhase(t1);  // Check for phase transitions

      // Update bucket stats every 5 seconds
      if((t1 - lastHotspotUpdate) > 5.0) {
        UpdateBucketStatistics(t1);
        CalculateHotspotScores();
        UpdateTopHotspots();
        lastHotspotUpdate = t1;
      }
    }

    // Display stats
    if(isAlive(params) && !endOfSearch) {
      // Calculate tame/wild ratio (1.000 = 50/50)
      double twRatio = (wildCount > 0) ? ((double)tameCount / (double)wildCount) : 0.0;

      // Calculate compact gap values - convert to billions
      // Full 128-bit value: (high * 2^64 + low) / 1e9
      double gap128 = (double)lastGap.i64[1] * 18446744073709551616.0 + (double)lastGap.i64[0];
      double lowestGap128 = (double)lowestGap.i64[1] * 18446744073709551616.0 + (double)lowestGap.i64[0];
      double currentGap = gap128 / 1000000000.0;
      // Check if lowestGap is still at initial max value (not yet set)
      bool lowestGapValid = !(lowestGap.i64[0] == 0xFFFFFFFFFFFFFFFFULL &&
                              lowestGap.i64[1] == 0x3FFFFFFFFFFFFFFFULL);
      double lowest = lowestGapValid ? (lowestGap128 / 1000000000.0) : 0.0;

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

      printf("\n\033[K[LZB:%u/%.0f][P:%.1f%%][50%%:%s][90%%:%s]%s  ",
        localMaxLZB,
        expectedMaxLZB,
        collisionProb * 100.0,
        (eta50 >= 0) ? GetTimeStr(eta50).c_str() : "N/A",
        (eta90 >= 0) ? GetTimeStr(eta90).c_str() : "N/A",
        hotBucketsStr
        );

      // Third line - Graduated DP phase info
      if(gradConfig.enabled) {
        char phaseInfo[128];
        GetPhaseInfo(phaseInfo, sizeof(phaseInfo));
        uint32_t topHotspotCount = (uint32_t)topHotspots.size();
        uint32_t topBucketId = topHotspots.empty() ? 0 : topHotspots[0].bucketId;
        double topScore = topHotspots.empty() ? 0.0 : topHotspots[0].score;

        printf("\n\033[K[PHASE:%s][Hotspots:%u][Top:0x%X(%.1f)]  \033[F\033[F",
          phaseInfo,
          topHotspotCount,
          topBucketId,
          topScore
          );
      } else {
        // Just move cursor back up if graduated DP is disabled
        printf("\033[F");
      }

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

