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

#include "HashTable.h"
#include <stdio.h>
#include <math.h>
#ifndef WIN64
#include <string.h>
#endif

HashTable::HashTable() {

  memset(E,0,sizeof(E));
  
}

void HashTable::Reset() {

  for(uint32_t h = 0; h < HASH_SIZE; h++) {
    safe_free(E[h].items);
    safe_free(E[h].states);
    E[h].maxItem = 0;
    E[h].nbItem = 0;
    E[h].nbTombstone = 0;
  }

}

uint64_t HashTable::GetNbItem() {

  uint64_t totalItem = 0;
  for(uint64_t h = 0; h < HASH_SIZE; h++) 
    totalItem += (uint64_t)E[h].nbItem;

  return totalItem;

}

void HashTable::Allocate(uint64_t h,uint32_t size) {

  E[h].maxItem = size;
  E[h].nbItem = 0;
  E[h].nbTombstone = 0;
  E[h].items = (ENTRY *)malloc(sizeof(ENTRY) * size);
  E[h].states = (uint8_t *)malloc(sizeof(uint8_t) * size);
  memset(E[h].states,SLOT_EMPTY,size);

}

void HashTable::Convert(Int *x,Int *d,uint32_t type,uint64_t *h,int128_t *X,int128_t *D) {

  uint64_t sign = 0;
  uint64_t type64 = (uint64_t)type << 62;

  X->i64[0] = x->bits64[0];
  X->i64[1] = x->bits64[1];

  // Probability of failure (1/2^128)
  if(d->bits64[3] > 0x7FFFFFFFFFFFFFFFULL) {
    Int N(d);
    N.ModNegK1order();
    D->i64[0] = N.bits64[0];
    D->i64[1] = N.bits64[1] & 0x3FFFFFFFFFFFFFFFULL;
    sign = 1ULL << 63;
  } else {
    D->i64[0] = d->bits64[0];
    D->i64[1] = d->bits64[1] & 0x3FFFFFFFFFFFFFFFULL;
  }

  D->i64[1] |= sign;
  D->i64[1] |= type64;

  *h = (x->bits64[2] & HASH_MASK);

}


#define AV1() if(pnb1) { ::fread(&e1,32,1,f1); pnb1--; }
#define AV2() if(pnb2) { ::fread(&e2,32,1,f2); pnb2--; }

int HashTable::MergeH(uint32_t h,FILE* f1,FILE* f2,FILE* fd,uint32_t* nbDP,uint32_t *duplicate,Int* d1,uint32_t* k1,Int* d2,uint32_t* k2) {

  // Merge by line
  // N comparison but avoid slow item allocation
  // return ADD_OK or ADD_COLLISION if a COLLISION is detected

  uint32_t nb1;
  uint32_t m1;
  uint32_t nb2;
  uint32_t m2;
  *duplicate = 0;
  *nbDP = 0;

  ::fread(&nb1,sizeof(uint32_t),1,f1);
  ::fread(&m1,sizeof(uint32_t),1,f1);
  ::fread(&nb2,sizeof(uint32_t),1,f2);
  ::fread(&m2,sizeof(uint32_t),1,f2);

  // Maximum in destination
  uint32_t nbd = 0;
  uint32_t md = nb1 + nb2;

  if(md==0) {

    ::fwrite(&md,sizeof(uint32_t),1,fd);
    ::fwrite(&md,sizeof(uint32_t),1,fd);
    return ADD_OK;

  }

  ENTRY *output = (ENTRY *)malloc( md * sizeof(ENTRY) );

  ENTRY e1;
  ENTRY e2;

  uint32_t pnb1 = nb1;
  uint32_t pnb2 = nb2;
  AV1();
  AV2();
  bool end1 = (nb1 == 0);
  bool end2 = (nb2 == 0);
  bool collisionFound = false;

  while(!(end1 && end2)) {

    if( !end1 && !end2 ) {

      int comp = compare(&e1.x,&e2.x);
      if(comp < 0) {
        memcpy(output+nbd,&e1,32);
        nbd++;
        AV1();
        nb1--;
      } else if (comp==0) {
        if((e1.d.i64[0] == e2.d.i64[0]) && (e1.d.i64[1] == e2.d.i64[1])) {
          *duplicate = *duplicate + 1;
        } else {
          // Collision
          CalcDistAndType(e1.d,d1,k1);
          CalcDistAndType(e2.d,d2,k2);
          collisionFound = true;
        }
        memcpy(output + nbd,&e1,32);
        nbd++;
        AV1();
        AV2();
        nb1--;
        nb2--;
      } else {
        memcpy(output + nbd,&e2,32);
        nbd++;
        AV2();
        nb2--;
      }

    } else if( !end1 && end2 ) {

      memcpy(output + nbd,&e1,32);
      nbd++;
      AV1();
      nb1--;

    } else if( end1 && !end2) {

      memcpy(output + nbd,&e2,32);
      nbd++;
      AV2();
      nb2--;

    }

    end1 = (nb1 == 0);
    end2 = (nb2 == 0);

  }

  // write output

  // Round md to next multiple of 4
  if(nbd%4==0) {
    md = nbd;
  } else {
    md = ((nbd/4)+1)*4;
  }

  ::fwrite(&nbd,sizeof(uint32_t),1,fd);
  ::fwrite(&md,sizeof(uint32_t),1,fd);
  ::fwrite(output,32,nbd,fd);
  free(output);

  *nbDP = nbd;
  return (collisionFound?ADD_COLLISION:ADD_OK);

}

int HashTable::Add(Int *x,Int *d,uint32_t type) {

  int128_t X;
  int128_t D;
  uint64_t h;
  Convert(x,d,type,&h,&X,&D);
  ENTRY e;
  e.x = X;
  e.d = D;
  return Add(h,e);

}

bool HashTable::NeedRehash(const HASH_ENTRY& he) {

  // Keep load factor under 70%, taking tombstones into account
  uint64_t used = (uint64_t)he.nbItem + (uint64_t)he.nbTombstone + 1ULL;
  return used * 10 > (uint64_t)he.maxItem * 7;

}

int HashTable::Add(uint64_t h,int128_t *x,int128_t *d) {

  ENTRY e;
  e.x = *x;
  e.d = *d;
  return Add(h,e);

}

void HashTable::CalcDistAndType(int128_t d,Int* kDist,uint32_t* kType) {

  *kType = (d.i64[1] & 0x4000000000000000ULL) != 0;
  int sign = (d.i64[1] & 0x8000000000000000ULL) != 0;
  d.i64[1] &= 0x3FFFFFFFFFFFFFFFULL;

  kDist->SetInt32(0);
  kDist->bits64[0] = d.i64[0];
  kDist->bits64[1] = d.i64[1];
  if(sign) kDist->ModNegK1order();

}

uint32_t HashTable::StartIndex(const ENTRY& e,const HASH_ENTRY& he) {

  uint64_t key = e.x.i64[0] ^ e.x.i64[1];
  return (uint32_t)(key % he.maxItem);

}

uint32_t HashTable::InsertInBucket(HASH_ENTRY& he,const ENTRY& e,bool detectCollision) {

  uint32_t start = StartIndex(e,he);
  uint32_t insertPos = he.maxItem;

  for(uint32_t p = 0; p < he.maxItem; p++) {
    uint32_t pos = (start + p);
    if(pos >= he.maxItem) pos -= he.maxItem;

    uint8_t state = he.states[pos];
    if(state == SLOT_EMPTY) {
      if(insertPos == he.maxItem) insertPos = pos;
      he.items[insertPos] = e;
      he.states[insertPos] = SLOT_OCCUPIED;
      he.nbItem++;
      if(insertPos != pos) he.nbTombstone--;
      return ADD_OK;
    }

    if(state == SLOT_TOMBSTONE) {
      if(insertPos == he.maxItem) insertPos = pos;
      continue;
    }

    if(detectCollision) {
      int comp = compare((int128_t *)&e.x,(int128_t *)&he.items[pos].x);
      if(comp == 0) {
        if((e.d.i64[0] == he.items[pos].d.i64[0]) && (e.d.i64[1] == he.items[pos].d.i64[1])) {
          // Same point added 2 times or collision in same herd !
          return ADD_DUPLICATE;
        }

        // Collision
        CalcDistAndType(he.items[pos].d , &kDist, &kType);
        return ADD_COLLISION;
      }
    }
  }

  return ADD_COLLISION;

}

void HashTable::ReHash(uint64_t h,uint32_t newSize) {

  HASH_ENTRY old = E[h];
  Allocate(h,newSize);

  for(uint32_t i = 0; i < old.maxItem; i++) {
    if(old.states && old.states[i] == SLOT_OCCUPIED) {
      InsertInBucket(E[h],old.items[i],false);
    }
  }

  safe_free(old.items);
  safe_free(old.states);

}

int HashTable::Add(uint64_t h,const ENTRY& e) {

  if(E[h].maxItem == 0) {
    Allocate(h,16);
  }

  if(NeedRehash(E[h])) {
    ReHash(h,E[h].maxItem * 2);
  }

  int res = InsertInBucket(E[h],e,true);
  if(res == ADD_COLLISION && NeedRehash(E[h])) {
    ReHash(h,E[h].maxItem * 2);
    res = InsertInBucket(E[h],e,true);
  }

  return res;

}

int HashTable::compare(int128_t *i1,int128_t *i2) {

  uint64_t *a = i1->i64;
  uint64_t *b = i2->i64;

  if(a[1] == b[1]) {
    if(a[0] == b[0]) {
      return 0;
    } else {
      return (a[0] > b[0]) ? 1 : -1;
    }
  } else {
    return (a[1] > b[1]) ? 1 : -1;
  }

}

std::string HashTable::GetSizeInfo() {

  char *unit;
  uint64_t totalByte = sizeof(E);
  uint64_t usedByte = HASH_SIZE*3*sizeof(uint32_t);

  for (int h = 0; h < HASH_SIZE; h++) {
    totalByte += sizeof(ENTRY) * E[h].maxItem;
    totalByte += sizeof(uint8_t) * E[h].maxItem;
    usedByte += sizeof(ENTRY) * E[h].nbItem;
    usedByte += sizeof(uint8_t) * E[h].maxItem;
  }

  unit = "MB";
  double totalMB = (double)totalByte / (1024.0*1024.0);
  double usedMB = (double)usedByte / (1024.0*1024.0);
  if(totalMB > 1024) {
    totalMB /= 1024;
    usedMB /= 1024;
    unit = "GB";
  }
  if(totalMB > 1024) {
    totalMB /= 1024;
    usedMB /= 1024;
    unit = "TB";
  }

  char ret[256];
  sprintf(ret,"%.1f/%.1f%s",usedMB,totalMB,unit);

  return std::string(ret);

}

std::string HashTable::GetStr(int128_t *i) {

  std::string ret;
  char tmp[256];
  for(int n=3;n>=0;n--) {
    ::sprintf(tmp,"%08X",i->i32[n]); 
    ret += std::string(tmp);
  }
  return ret;

}

void HashTable::SaveTable(FILE* f) {
  SaveTable(f,0,HASH_SIZE,true);
}

void HashTable::SaveTable(FILE* f,uint32_t from,uint32_t to,bool printPoint) {

  uint64_t point = GetNbItem() / 16;
  uint64_t pointPrint = 0;

  for(uint32_t h = from; h < to; h++) {
    fwrite(&E[h].nbItem,sizeof(uint32_t),1,f);
    fwrite(&E[h].maxItem,sizeof(uint32_t),1,f);
    if(E[h].states) {
      for(uint32_t i = 0; i < E[h].maxItem; i++) {
        if(E[h].states[i] == SLOT_OCCUPIED) {
          fwrite(&(E[h].items[i].x),16,1,f);
          fwrite(&(E[h].items[i].d),16,1,f);
          if(printPoint) {
            pointPrint++;
            if(pointPrint > point) {
              ::printf(".");
              pointPrint = 0;
            }
          }
        }
      }
    }
  }

}

void HashTable::SeekNbItem(FILE* f,bool restorePos) {

  Reset();

#ifdef WIN64
  uint64_t org = (uint64_t)_ftelli64(f);
#else
  uint64_t org = (uint64_t)ftello(f);
#endif

  SeekNbItem(f,0,HASH_SIZE);

  if( restorePos ) {
    // Restore position
#ifdef WIN64
    _fseeki64(f,org,SEEK_SET);
#else
    fseeko(f,org,SEEK_SET);
#endif
  }

}

void HashTable::SeekNbItem(FILE* f,uint32_t from,uint32_t to) {

  for(uint32_t h = from; h < to; h++) {

    fread(&E[h].nbItem,sizeof(uint32_t),1,f);
    fread(&E[h].maxItem,sizeof(uint32_t),1,f);

    uint64_t hSize = 32ULL * E[h].nbItem;
#ifdef WIN64
    _fseeki64(f,hSize,SEEK_CUR);
#else
    fseeko(f,hSize,SEEK_CUR);
#endif

  }

}

void HashTable::LoadTable(FILE* f,uint32_t from,uint32_t to) {

  Reset();

  for(uint32_t h = from; h < to; h++) {

    uint32_t nbItem;
    fread(&nbItem,sizeof(uint32_t),1,f);
    fread(&E[h].maxItem,sizeof(uint32_t),1,f);

    if(E[h].maxItem > 0)
      Allocate(h,E[h].maxItem);

    for(uint32_t i = 0; i < nbItem; i++) {
      ENTRY e;
      fread(&(e.x),16,1,f);
      fread(&(e.d),16,1,f);
      InsertInBucket(E[h],e,false);
    }

  }


}

void HashTable::LoadTable(FILE *f) {

  LoadTable(f,0,HASH_SIZE);

}

void HashTable::PrintInfo() {

  uint16_t max = 0;
  uint32_t maxH = 0;
  uint16_t min = 65535;
  uint32_t minH = 0;
  double std = 0;
  double avg = (double)GetNbItem() / (double)HASH_SIZE;

  for(uint32_t h=0;h<HASH_SIZE;h++) {
    if(E[h].nbItem>max) {
      max= E[h].nbItem;
      maxH = h;
    }
    if(E[h].nbItem<min) {
      min= E[h].nbItem;
      minH = h;
    }
    std += (avg - (double)E[h].nbItem)*(avg - (double)E[h].nbItem);
  }
  std /= (double)HASH_SIZE;
  std = sqrt(std);

  uint64_t count = GetNbItem();

  ::printf("DP Size   : %s\n",GetSizeInfo().c_str());
#ifdef WIN64
  ::printf("DP Count  : %I64d 2^%.3f\n",count,log2((double)count));
#else
  ::printf("DP Count  : %" PRId64 " 2^%.3f\n",count,log2(count));
#endif
  ::printf("HT Max    : %d [@ %06X]\n",max,maxH);
  ::printf("HT Min    : %d [@ %06X]\n",min,minH);
  ::printf("HT Avg    : %.2f \n",avg);
  ::printf("HT SDev   : %.2f \n",std);

  //for(int i=0;i<(int)E[maxH].nbItem;i++) {
  //  ::printf("[%2d] %s\n",i,GetStr(&E[maxH].items[i]->x).c_str());
  //  ::printf("[%2d] %s\n",i,GetStr(&E[maxH].items[i]->d).c_str());
  //}

}
