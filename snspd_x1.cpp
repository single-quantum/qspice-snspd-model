//    cl /LD /std:c++17 snspd_x1.cpp kernel32.lib

#include <malloc.h>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <execution>
#include <vector>

const double TSUBERROR = 1.01;
const double GAMMA = 240;
const double ALPHA_SF = 800;
const double LORENTZ = 2.45e-8;
const double G = 9.8;

union uData
{
   bool b;
   char c;
   unsigned char uc;
   short s;
   unsigned short us;
   int i;
   unsigned int ui;
   float f;
   double d;
   long long int i64;
   unsigned long long int ui64;
   char *str;
   unsigned char *bytes;
};

struct sSNSPD_X1
{
    double resistance;
    //double resistancePrev;
    double* temperatures;
    double  width;
    double  length;
    double  thickness;
    double  dlength;
    int     resolution;
    double  Tc;
    double  Tsub;
    double  Ic0K;
    double  Rsheet;
    double  Rsegment;
    int     photonnumber;
    bool    hotspot;
    double  ths;
    double  Ths;
    double  sizehs;
    double  time;
    double* diagonal;
    double* off_diagonal;
    double* right_hand_side;
    double  Tmax;
    bool    running;
    double prevCurrent;
    double Lkin;
    double V1;
    double V2;
    double prevR;
    double starttime;
    double Rmin;
    int mod;
    double A;
    double y;

    // Constructor (optional, but good practice in C++)
    sSNSPD_X1(union uData *data) : resistance(data[16].d), temperatures(nullptr), time(0),
    width(data[2].d), length(data[3].d), thickness(data[4].d), resolution(data[5].i), Tsub(data[7].d),
    Tc(data[6].d),  Ic0K(data[8].d), Rsheet(data[9].d), hotspot(data[10].b), ths(data[11].d), 
    Ths(data[12].d), sizehs(data[13].d), photonnumber(data[14].i), Lkin(data[15].d), Rmin(data[16].d), mod(data[17].i){
        temperatures = new double[resolution]; // Use '->' to access member via pointer
        diagonal = new double[resolution];
        off_diagonal = new double[resolution];
        right_hand_side = new double[resolution];

        A = 2.43*Tc*GAMMA;

        std::fill(temperatures, temperatures + resolution, Tsub);
        dlength    = length / (double)resolution;
        Rsegment   = Rsheet * dlength / width;
        Tmax       = Tsub;
        running    = false;

    }

    // Destructor (essential for managing dynamically allocated memory)
    ~sSNSPD_X1() {
        delete[] temperatures;
        delete[] diagonal;
        delete[] off_diagonal;
        delete[] right_hand_side;
    }
};

// int DllMain() must exist and return 1 for a process to load the .DLL
// See https://docs.microsoft.com/en-us/windows/win32/dlls/dllmain for more information.
int __stdcall DllMain(void *module, unsigned int reason, void *reserved) { return 1; }

// #undef pin names lest they collide with names in any header file(s) you might include.
#undef IN1
#undef IN2
#undef ROUT



// All Thermal functions
double calcAlpha(double temperature){
    return ALPHA_SF * pow(temperature,3);
}

double calcKapa(double thickness, double Rsheet, double temperature){
    return ( LORENTZ * temperature ) / (Rsheet * thickness);
}

double calcKapaSC(double thickness, double Rsheet, double Tc, double temperature){
    return calcKapa(thickness, Rsheet, temperature) * temperature / Tc;
}

double calcPhononHC(double temperature){
    return G * pow(temperature,3);
}

double calcPow( double Tc,  double T){
    return -2.15*Tc*(1-pow(T/Tc,2))/T;
}

double calcElectronHCSC(double Tc, double temperature, double A){
    return  A * exp(calcPow(Tc, temperature));
}

double calcElectronHC(double temperature){
    return temperature * GAMMA;
}

double calcIc(double Ic0K, double Tc, double temperature){
    return Ic0K * pow(1.0-pow(temperature / Tc, 2),2);
}

bool isSC(double current, double temperature, double Ic0K , double Tc){
    return fabs(current) < calcIc(Ic0K, Tc, temperature) && temperature <= Tc;
}

void createHotspot(sSNSPD_X1 *opaque, double current){
    int PNR = opaque->photonnumber;
    int hotspot_segments = (int)(opaque->sizehs * opaque->resolution / opaque->length);
    double resistance = opaque->Rmin;
    double Tmax = opaque->Tsub;

    for (int n = 0; n < PNR; ++n) {
        int start_index_hotspot = (1 + n)*opaque->resolution / (1 + PNR) - (hotspot_segments / 2);
        for (int i = start_index_hotspot; i < start_index_hotspot + hotspot_segments; ++i) {
            opaque->temperatures[i] = opaque->Ths;
            if (!isSC(current, opaque->temperatures[i], opaque->Ic0K, opaque->Tc)){
                resistance = opaque->Rsegment ++;
            }
        }
    }
    opaque->resistance = resistance;
    opaque->Tmax = opaque->Ths;

}

void diagonalSC(sSNSPD_X1 *opaque, int index, double dt, double* diagonal, double* off_diagonal, double* right_hand_side){
    double alpha = calcAlpha(opaque->temperatures[index]);
    double kapa = calcKapaSC(opaque->thickness, opaque->Rsheet, opaque->Tc, opaque->temperatures[index]);
    double heat_cap = calcElectronHCSC(opaque->Tc, opaque->temperatures[index], opaque->A) + calcPhononHC(opaque->temperatures[index]);

    double r = kapa * dt / (2 * pow(opaque->dlength,2) * heat_cap);
    double h = alpha * dt / (2 * opaque->thickness * heat_cap);
    double g = dt / heat_cap * ( opaque->Tsub * alpha / opaque->thickness);

    off_diagonal[index] = -r;
    diagonal[index]=  1 + h + 2 * r;
    right_hand_side[index] = opaque->temperatures[index] * (1 - h - 2 * r) + r * (opaque->temperatures[index+1] + opaque->temperatures[index-1]) + g;
}

void  diagonalNSC(sSNSPD_X1 *opaque, int index, double dt, double current, double* diagonal, double* off_diagonal, double* right_hand_side){
    double alpha = calcAlpha(opaque->temperatures[index]);
    double kapa = calcKapa(opaque->thickness, opaque->Rsheet, opaque->temperatures[index]);
    double heat_cap = calcElectronHC(opaque->temperatures[index]) + calcPhononHC(opaque->temperatures[index]);

    double r = kapa * dt / (2 * pow(opaque->dlength,2) * heat_cap);
    double h = alpha * dt / (2 * opaque->thickness * heat_cap);
    double g = dt / heat_cap * ( opaque->Tsub * alpha  + pow(current / opaque->width, 2) * opaque->Rsheet) / opaque->thickness;

    off_diagonal[index] = -r;
    diagonal[index]=  1 + h + 2 * r;
    right_hand_side[index] = opaque->temperatures[index] * (1 - h - 2 * r) + r * (opaque->temperatures[index+1] + opaque->temperatures[index-1]) + g;
}

void calcTotalResitance(sSNSPD_X1 *opaque, double current, double dt){
    double resistance = 0;
    double Tmax = opaque->Tsub;

    opaque->diagonal[0] = opaque->diagonal[opaque->resolution - 1] = 1.0;
    opaque->right_hand_side[0] = opaque->right_hand_side[opaque->resolution - 1] = opaque->Tsub;
    opaque->off_diagonal[0] = opaque->off_diagonal[opaque->resolution - 1] = 0;

    std::vector<int> indices(opaque->resolution - 2);
    std::iota(indices.begin(), indices.end(), 1);

    std::for_each(std::execution::par_unseq, indices.begin(), indices.end(),
                  [&](int i) {
                      if (isSC(current, opaque->temperatures[i], opaque->Ic0K, opaque->Tc)){
                          diagonalSC(opaque, i, dt, opaque->diagonal, opaque->off_diagonal, opaque->right_hand_side);
                      } else {
                          diagonalNSC(opaque, i, dt, current, opaque->diagonal, opaque->off_diagonal, opaque->right_hand_side);
                      }
                  });
    double m;
    double ai;
    double ci1;
    double* x;
    x = new double[opaque->resolution];

    for (int it = 1; it < opaque->resolution; ++it) {
        ai = opaque->off_diagonal[it];
        ci1 = opaque->off_diagonal[it-1];
        m = ai / opaque->diagonal[it-1];
        opaque->diagonal[it] = opaque->diagonal[it] - m * ci1;
        opaque->right_hand_side[it] = opaque->right_hand_side[it] - m * opaque->right_hand_side[it-1];
    }
    x[opaque->resolution - 1] = opaque->right_hand_side[opaque->resolution - 1] / opaque->diagonal[opaque->resolution -1];
    for (int il = opaque->resolution-2; il >=1; --il){
        x[il] = (opaque->right_hand_side[il] - opaque->off_diagonal[il] *  x[il + 1]) / opaque->diagonal[il];
        if (!isSC(current, x[il], opaque->Ic0K, opaque->Tc)){
            resistance += opaque->Rsegment;
        }
    }
    opaque->temperatures =x;

    opaque->resistance = resistance;

}

double calcCurrentTrapezoidalRule(double V1, double V2, double V1prev, double V2prev, double Lkin, double dt, double R, double Rprev, double current){
   double dV = (V1-V2);
   double dVprev = (V1prev-V2prev);
   double RL = 2.0*Lkin/dt;
   return ( (RL-Rprev)* current + dV + dVprev)/(R+RL);
}

double calcBackwardsEuler(double V1, double V2,double Lkin, double dt, double R,double current){
   double dV = (V1-V2);
   double RL = Lkin/dt;
   return (RL* current + dV)/(R+RL);
}

extern "C" __declspec(dllexport) void snspd_x1(struct sSNSPD_X1 **opaque, double t, union uData *data)
{
   double IN1          = data[ 0].d; // input
   double IN2          = data[ 1].d; // input
   double &ROUT = data[18].d; // output
   //double &ROUTPREV = data[15].d; // output
   sSNSPD_X1 *inst = *opaque;

   if(!inst){
        inst = *opaque = new sSNSPD_X1(data);
        inst->starttime = t;
   }

   if (t==inst->starttime) {
        ROUT = 0;
        inst->time = t;
        return;
   }

   double dt = t - inst->time;
   inst->time = t;
   if (inst->hotspot && t >= inst->ths){
        inst->hotspot = false;
        createHotspot(inst, IN1);
   }
   else{
    calcTotalResitance(inst, IN1, dt);
   }

   ROUT = inst->resistance;
}

extern "C" __declspec(dllexport) void Destroy(struct sSNSPD_X1 *inst)
{
    delete inst;
    inst = nullptr;
}
