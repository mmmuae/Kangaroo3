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
#include "SECPK1/SECP256k1.h"
#include "SECPK1/Int.h"
#include "GPU/GPUEngine.h"
#include <fstream>
#include <string>
#include <string.h>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <vector>
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

using namespace std;

#define CHECKARG(opt,n) if(a>=argc-1) {::printf(opt " missing argument #%d\n",n);exit(0);} else {a++;}

// ------------------------------------------------------------------------------------------

void printUsage() {

  printf("Kangaroo [-v] [-t nbThread] [-d dpBit] [gpu] [-check]\n");
  printf("         [-gpuId gpuId1[,gpuId2,...]] [-g g1x,g1y[,g2x,g2y,...]]\n");
  printf("         inFile\n");
  printf(" -v: Print version\n");
  printf(" -gpu: Enable gpu calculation\n");
  printf(" -gpuId gpuId1,gpuId2,...: List of GPU(s) to use, default is 0\n");
  printf(" -g g1x,g1y,g2x,g2y,...: Specify GPU(s) kernel gridsize, default is 2*(MP),2*(Core/MP)\n");
  printf(" -d: Specify number of leading zeros for the DP method (default is auto)\n");
  printf(" -t nbThread: Secify number of thread\n");
  printf(" -w workfile: Specify file to save work into (current processed key only)\n");
  printf(" -i workfile: Specify file to load work from (current processed key only)\n");
  printf(" --start-dec <start_dec>  --end-dec <end_dec>  --pubkey <hex>\n");
  printf(" --start-hex <start_hex>  --end-hex <end_hex>  --pubkey <hex>\n");
  printf(" -wi workInterval: Periodic interval (in seconds) for saving work\n");
  printf(" -ws: Save kangaroos in the work file\n");
  printf(" -wss: Save kangaroos via the server\n");
  printf(" -wsplit: Split work file of server and reset hashtable\n");
  printf(" -wm file1 file2 destfile: Merge work file\n");
  printf(" -wmdir dir destfile: Merge directory of work files\n");
  printf(" -wt timeout: Save work timeout in millisec (default is 3000ms)\n");
  printf(" -winfo file1: Work file info file\n");
  printf(" -wpartcreate name: Create empty partitioned work file (name is a directory)\n");
  printf(" -wcheck worfile: Check workfile integrity\n");
  printf(" -m maxStep: number of operations before give up the search (maxStep*expected operation)\n");
  printf(" -s: Start in server mode\n");
  printf(" -c server_ip: Start in client mode and connect to server server_ip\n");
  printf(" -sp port: Server port, default is 17403\n");
  printf(" -nt timeout: Network timeout in millisec (default is 3000ms)\n");
  printf(" -o fileName: output result to fileName\n");
  printf(" -l: List cuda enabled devices\n");
  printf(" -check: Check GPU kernel vs CPU\n");
  printf(" inFile: intput configuration file\n");
  exit(0);

}

// ------------------------------------------------------------------------------------------

int getInt(string name,char *v) {

  int r;

  try {

    r = std::stoi(string(v));

  } catch(std::invalid_argument&) {

    printf("Invalid %s argument, number expected\n",name.c_str());
    exit(-1);

  }

  return r;

}

double getDouble(string name,char *v) {

  double r;

  try {

    r = std::stod(string(v));

  } catch(std::invalid_argument&) {

    printf("Invalid %s argument, number expected\n",name.c_str());
    exit(-1);

  }

  return r;

}

// ------------------------------------------------------------------------------------------

void getInts(string name,vector<int> &tokens,const string &text,char sep) {

  size_t start = 0,end = 0;
  tokens.clear();
  int item;

  try {

    while((end = text.find(sep,start)) != string::npos) {
      item = std::stoi(text.substr(start,end - start));
      tokens.push_back(item);
      start = end + 1;
    }

    item = std::stoi(text.substr(start));
    tokens.push_back(item);

  }
  catch(std::invalid_argument &) {

    printf("Invalid %s argument, number expected\n",name.c_str());
    exit(-1);

  }

}
// ------------------------------------------------------------------------------------------

// Default params
static int dp = -1;
static int nbCPUThread;
static string configFile = "";
static bool checkFlag = false;
static bool gpuEnable = false;
static vector<int> gpuId = { 0 };
static vector<int> gridSize;
static string workFile = "";
static string checkWorkFile = "";
static string iWorkFile = "";
static uint32_t savePeriod = 60;
static bool saveKangaroo = false;
static bool saveKangarooByServer = false;
static string merge1 = "";
static string merge2 = "";
static string mergeDest = "";
static string mergeDir = "";
static string infoFile = "";
static double maxStep = 0.0;
static int wtimeout = 3000;
static int ntimeout = 3000;
static int port = 17403;
static bool serverMode = false;
static string serverIP = "";
static string outputFile = "";
static bool splitWorkFile = false;

static bool isHexString(const string &value) {
  if(value.empty()) {
    return false;
  }
  for(char c : value) {
    if(!isxdigit(static_cast<unsigned char>(c))) {
      return false;
    }
  }
  return true;
}

static string stripHexPrefix(const string &value) {
  if(value.size() >= 2 && value[0] == '0' && (value[1] == 'x' || value[1] == 'X')) {
    return value.substr(2);
  }
  return value;
}

static void setIntBase10(Int &target, const string &value) {
  vector<char> buff(value.begin(), value.end());
  buff.push_back('\0');
  target.SetBase10(buff.data());
}

static void setIntBase16(Int &target, const string &value) {
  vector<char> buff(value.begin(), value.end());
  buff.push_back('\0');
  target.SetBase16(buff.data());
}

static std::string toHex64(Int value) {
  std::string hex = value.GetBase16();  // OK: value is non-const here
  if (hex.size() > 64) {
    std::printf("Range value exceeds 256-bit limit\n");
    std::exit(-1);
  }
  if (hex.size() < 64) {
    hex.insert(hex.begin(), static_cast<std::string::size_type>(64 - hex.size()), '0');
  }
  return hex;
}

class ScopedTempFile {
public:
  ScopedTempFile() : file(nullptr) {}
  ~ScopedTempFile() {
    cleanup();
  }

  void create() {
    cleanup();
#ifdef _WIN32
    char tempPath[MAX_PATH];
    if(GetTempPathA(MAX_PATH, tempPath) == 0) {
      throw std::runtime_error("GetTempPathA failed");
    }
    char tempFile[MAX_PATH];
    if(GetTempFileNameA(tempPath, "kang", 0, tempFile) == 0) {
      throw std::runtime_error("GetTempFileNameA failed");
    }
    path.assign(tempFile);
    file = fopen(path.c_str(), "w+");
    if(!file) {
      throw std::runtime_error("Unable to open temporary file");
    }
#else
    const char *tmpdir = getenv("TMPDIR");
    string base = tmpdir ? string(tmpdir) : string(P_tmpdir ? P_tmpdir : "/tmp");
    if(base.empty()) {
      base = "/tmp";
    }
    if(base.back() != '/') {
      base.push_back('/');
    }
    string templ = base + "kangarooXXXXXX";
    vector<char> buf(templ.begin(), templ.end());
    buf.push_back('\0');
    int fd = mkstemp(buf.data());
    if(fd == -1) {
      throw std::runtime_error("Unable to create temporary file");
    }
    path.assign(buf.data());
    file = fdopen(fd, "w+");
    if(!file) {
      close(fd);
      throw std::runtime_error("Unable to open temporary file stream");
    }
#endif
  }

  FILE *get() const {
    return file;
  }

  const string &getPath() const {
    return path;
  }

  void cleanup() {
    if(file) {
      fclose(file);
      file = nullptr;
    }
    if(!path.empty()) {
      std::remove(path.c_str());
      path.clear();
    }
  }

private:
  FILE *file;
  string path;
};

int main(int argc, char* argv[]) {

#ifdef USE_SYMMETRY
  printf("Kangaroo v" RELEASE " (with symmetry)\n");
#else
  printf("Kangaroo v" RELEASE "\n");
#endif

  // Global Init
  Timer::Init();
  rseed(Timer::getSeed32());

  // Init SecpK1
  Secp256K1 *secp = new Secp256K1();
  secp->Init();

  int a = 1;
  nbCPUThread = Timer::getCoreNumber();

  string startDecArg;
  string endDecArg;
  string startHexArg;
  string endHexArg;
  string pubkeyArg;

  while (a < argc) {

    if(strcmp(argv[a], "-t") == 0) {
      CHECKARG("-t",1);
      nbCPUThread = getInt("nbCPUThread",argv[a]);
      a++;
    } else if(strcmp(argv[a],"-d") == 0) {
      CHECKARG("-d",1);
      dp = getInt("dpSize",argv[a]);
      a++;
    } else if (strcmp(argv[a], "-h") == 0) {
      printUsage();
    } else if(strcmp(argv[a],"-l") == 0) {

#ifdef WITHGPU
      GPUEngine::PrintCudaInfo();
#else
      printf("GPU code not compiled, use -DWITHGPU when compiling.\n");
#endif
      exit(0);

    } else if(strcmp(argv[a],"-w") == 0) {
      CHECKARG("-w",1);
      workFile = string(argv[a]);
      a++;
    } else if(strcmp(argv[a],"-i") == 0) {
      CHECKARG("-i",1);
      iWorkFile = string(argv[a]);
      a++;
    } else if(strcmp(argv[a],"-wm") == 0) {
      CHECKARG("-wm",1);
      merge1 = string(argv[a]);
      CHECKARG("-wm",2);
      merge2 = string(argv[a]);
      a++;
      if(a<argc) {
        // classic merge
        mergeDest = string(argv[a]);
        a++;
      }
    } else if(strcmp(argv[a],"-wmdir") == 0) {
      CHECKARG("-wmdir",1);
      mergeDir = string(argv[a]);
      CHECKARG("-wmdir",2);
      mergeDest = string(argv[a]);
      a++;
    }  else if(strcmp(argv[a],"-wcheck") == 0) {
      CHECKARG("-wcheck",1);
      checkWorkFile = string(argv[a]);
      a++;
    }  else if(strcmp(argv[a],"-winfo") == 0) {
      CHECKARG("-winfo",1);
      infoFile = string(argv[a]);
      a++;
    } else if(strcmp(argv[a],"-o") == 0) {
      CHECKARG("-o",1);
      outputFile = string(argv[a]);
      a++;
    } else if(strcmp(argv[a],"-wi") == 0) {
      CHECKARG("-wi",1);
      savePeriod = getInt("savePeriod",argv[a]);
      a++;
    } else if(strcmp(argv[a],"-wt") == 0) {
      CHECKARG("-wt",1);
      wtimeout = getInt("timeout",argv[a]);
      a++;
    } else if(strcmp(argv[a],"-nt") == 0) {
      CHECKARG("-nt",1);
      ntimeout = getInt("timeout",argv[a]);
      a++;
    } else if(strcmp(argv[a],"-m") == 0) {
      CHECKARG("-m",1);
      maxStep = getDouble("maxStep",argv[a]);
      a++;
    } else if(strcmp(argv[a],"-ws") == 0) {
      a++;
      saveKangaroo = true;
    } else if(strcmp(argv[a],"-wss") == 0) {
      a++;
      saveKangarooByServer = true;
    } else if(strcmp(argv[a],"-wsplit") == 0) {
      a++;
      splitWorkFile = true;
    } else if(strcmp(argv[a],"-wpartcreate") == 0) {
      CHECKARG("-wpartcreate",1);
      workFile = string(argv[a]);
      Kangaroo::CreateEmptyPartWork(workFile);
      exit(0);
    } else if(strcmp(argv[a],"-s") == 0) {
      a++;
      serverMode = true;
    } else if(strcmp(argv[a],"-c") == 0) {
      CHECKARG("-c",1);
      serverIP = string(argv[a]);
      a++;
    } else if(strcmp(argv[a],"-sp") == 0) {
      CHECKARG("-sp",1);
      port = getInt("serverPort",argv[a]);
      a++;
    } else if(strcmp(argv[a],"-gpu") == 0) {
      gpuEnable = true;
      a++;
    } else if(strcmp(argv[a],"-gpuId") == 0) {
      CHECKARG("-gpuId",1);
      getInts("gpuId",gpuId,string(argv[a]),',');
      a++;
    } else if(strcmp(argv[a],"-g") == 0) {
      CHECKARG("-g",1);
      getInts("gridSize",gridSize,string(argv[a]),',');
      a++;
    } else if(strcmp(argv[a],"-v") == 0) {
      ::exit(0);
    } else if(strcmp(argv[a],"-check") == 0) {
      checkFlag = true;
      a++;
    } else if(strcmp(argv[a],"--start-dec") == 0) {
      CHECKARG("--start-dec",1);
      startDecArg = string(argv[a]);
      a++;
    } else if(strcmp(argv[a],"--end-dec") == 0) {
      CHECKARG("--end-dec",1);
      endDecArg = string(argv[a]);
      a++;
    } else if(strcmp(argv[a],"--start-hex") == 0) {
      CHECKARG("--start-hex",1);
      startHexArg = stripHexPrefix(string(argv[a]));
      a++;
    } else if(strcmp(argv[a],"--end-hex") == 0) {
      CHECKARG("--end-hex",1);
      endHexArg = stripHexPrefix(string(argv[a]));
      a++;
    } else if(strcmp(argv[a],"--pubkey") == 0) {
      CHECKARG("--pubkey",1);
      pubkeyArg = stripHexPrefix(string(argv[a]));
      a++;
    } else if(a == argc - 1) {
      configFile = string(argv[a]);
      a++;
    } else {
      printf("Unexpected %s argument\n",argv[a]);
      exit(-1);
    }

  }

  bool hasDecRange = (!startDecArg.empty() || !endDecArg.empty());
  bool hasHexRange = (!startHexArg.empty() || !endHexArg.empty());
  bool hasRangeFlags = hasDecRange || hasHexRange;

  if(hasDecRange && hasHexRange) {
    printf("Specify either decimal or hex range flags, not both\n");
    exit(-1);
  }

  if(hasDecRange && (startDecArg.empty() || endDecArg.empty())) {
    printf("Both --start-dec and --end-dec must be provided\n");
    exit(-1);
  }

  if(hasHexRange && (startHexArg.empty() || endHexArg.empty())) {
    printf("Both --start-hex and --end-hex must be provided\n");
    exit(-1);
  }

  if(hasRangeFlags && pubkeyArg.empty()) {
    printf("--pubkey is required when range flags are used\n");
    exit(-1);
  }

  if(!pubkeyArg.empty() && !isHexString(pubkeyArg)) {
    printf("--pubkey must be hex encoded\n");
    exit(-1);
  }

  string normalizedStartHex;
  string normalizedEndHex;
  string normalizedPubkey;
  if(!pubkeyArg.empty()) {
    normalizedPubkey = pubkeyArg;
    std::transform(normalizedPubkey.begin(), normalizedPubkey.end(), normalizedPubkey.begin(), [](unsigned char c) {
      return static_cast<char>(toupper(c));
    });
  }

  if(hasRangeFlags) {
    Int startInt;
    Int endInt;

    if(hasDecRange) {
      for(char c : startDecArg) {
        if(!isdigit(static_cast<unsigned char>(c))) {
          printf("--start-dec must be a non-negative decimal\n");
          exit(-1);
        }
      }
      for(char c : endDecArg) {
        if(!isdigit(static_cast<unsigned char>(c))) {
          printf("--end-dec must be a non-negative decimal\n");
          exit(-1);
        }
      }
      setIntBase10(startInt, startDecArg);
      setIntBase10(endInt, endDecArg);
    } else {
      if(!isHexString(startHexArg) || !isHexString(endHexArg)) {
        printf("--start-hex and --end-hex must be hex encoded\n");
        exit(-1);
      }
      setIntBase16(startInt, startHexArg);
      setIntBase16(endInt, endHexArg);
    }

    if(startInt.IsGreater(&endInt)) {
      printf("Start range must be less than or equal to end range\n");
      exit(-1);
    }

    normalizedStartHex = toHex64(startInt);
    normalizedEndHex = toHex64(endInt);
  }

  if(gridSize.size() == 0) {
    for(int i = 0; i < gpuId.size(); i++) {
      gridSize.push_back(0);
      gridSize.push_back(0);
    }
  } else if(gridSize.size() != gpuId.size() * 2) {
    printf("Invalid gridSize or gpuId argument, must have coherent size\n");
    exit(-1);
  }

  Kangaroo *v = new Kangaroo(secp,dp,gpuEnable,workFile,iWorkFile,savePeriod,saveKangaroo,saveKangarooByServer,
                             maxStep,wtimeout,port,ntimeout,serverIP,outputFile,splitWorkFile);
  ScopedTempFile tempFile;

  if(checkFlag) {
    v->Check(gpuId,gridSize);
    exit(0);
  } else {
    if(checkWorkFile.length() > 0) {
      v->CheckWorkFile(nbCPUThread,checkWorkFile);
      exit(0);
    } if(infoFile.length()>0) {
      v->WorkInfo(infoFile);
      exit(0);
    } else if(mergeDir.length() > 0) {
      v->MergeDir(mergeDir,mergeDest);
      exit(0);
    } else if(merge1.length()>0) {
      v->MergeWork(merge1,merge2,mergeDest);
      exit(0);
    }

    if(hasRangeFlags) {
      try {
        tempFile.create();
      } catch(const std::exception &ex) {
        printf("Failed to create temporary configuration file: %s\n", ex.what());
        exit(-1);
      }

      FILE *fp = tempFile.get();
      if(!fp) {
        printf("Failed to open temporary configuration file\n");
        exit(-1);
      }

      if(fprintf(fp, "%s\n%s\n%s\n", normalizedStartHex.c_str(), normalizedEndHex.c_str(), normalizedPubkey.c_str()) < 0) {
        printf("Failed to write temporary configuration file\n");
        exit(-1);
      }
      fflush(fp);

      configFile = tempFile.getPath();
    }

    if(iWorkFile.length()>0) {
      if( !v->LoadWork(iWorkFile) ) {
        if(hasRangeFlags) {
          tempFile.cleanup();
        }
        exit(-1);
      }
    } else if(configFile.length()>0) {
      bool parsed = v->ParseConfigFile(configFile);
      if(hasRangeFlags) {
        tempFile.cleanup();
      }
      if( !parsed )
        exit(-1);
    } else {
      if(serverIP.length()==0) {
        ::printf("No input file to process\n");
        if(hasRangeFlags) {
          tempFile.cleanup();
        }
        exit(-1);
      }
    }
    if(serverMode)
      v->RunServer();
    else
      v->Run(nbCPUThread,gpuId,gridSize);
  }

  return 0;

}
