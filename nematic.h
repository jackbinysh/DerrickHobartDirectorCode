#include <string>
/* internal global variables */
extern int n,LL;
extern double *nx,*ny,*nz;
extern double *hx,*hy,*hz;
extern double K1mK2,K3mK2;
extern double q0;
extern std::string RunName;

// summary stats functions
void SplayTwistBendDensities(int j, double& splaysq, double &twist, double &bend, double &twistsq);
void writeStatistics(FILE* fileStats);  // output measurements


int i_vr(int l);
int j_vr(int l);
int k_vr(int l);
