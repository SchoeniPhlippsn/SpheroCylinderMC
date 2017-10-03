
//hardWalls main
const double rlist = 3.1;
const double rlist_2 = rlist*rlist;
const double maxpos = 0.2*rlist;
int sx,sy,sz;
long N; //# of particles
long seed;
long finstep;
long savestep;
long backupstep;
long acc;
double rho0; //density
double Vratio;
double Rdist;
double rdist;
double height;
double Rsp;
double Rsp_cube;
double rsp;
double rsp_cube;
double rsphere;
double pos_lambda;
double ori_lambda;
double vproc;
double dr, r;
std::vector<double> rij;
int randP;
double acceptance;

bool restart;
bool initCompress;

std::string file_name;
std::string RHO;
std::string NC;
std::string NS;
std::string VR;

#include "System.h"
//hardWalls main
class system Config;
class teilchen MovedParticle;
