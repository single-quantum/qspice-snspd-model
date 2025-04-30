// Automatically generated C++ file on Thu Apr 24 13:03:17 2025
//
// To build with Digital Mars C++ Compiler:
//
//    dmc -mn -WD snspd_x1.cpp kernel32.lib

#include <malloc.h>
#include <math.h>
#include <iostream>
#include <cstring> // For bzero (though not the C++ way)
#include <stdio.h>
//#include <print>
//#include <std>
//using namespace std;
const double MINRES = 1e-12;
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
    bool    hotspot;
    double  ths;
    double  Ths;
    double  sizehs;
    double  time;

    // Constructor (optional, but good practice in C++)
    sSNSPD_X1() : resistance(MINRES), temperatures(nullptr), time(0) {}

    // Destructor (essential for managing dynamically allocated memory)
    ~sSNSPD_X1() {
        delete[] temperatures;
    }
};

// int DllMain() must exist and return 1 for a process to load the .DLL
// See https://docs.microsoft.com/en-us/windows/win32/dlls/dllmain for more information.
int __stdcall DllMain(void *module, unsigned int reason, void *reserved) { return 1; }

void bzero(void *ptr, unsigned int count)
{
   unsigned char *first = (unsigned char *) ptr;
   unsigned char *last  = first + count;
   while(first < last)
      *first++ = '\0';
}

// #undef pin names lest they collide with names in any header file(s) you might include.
#undef IN
#undef OUT
const double HEATCAP_ELECTRON_SF = 240;
const double PHONON_HEATCAP_SF = 9.8;
const double ALPHA_SF = 800;
const double LORENTZ = 2.45e-8;


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
    return PHONON_HEATCAP_SF * pow(temperature,3);
}

double calcElectronHCSC(double Tc, double temperature){
    return 2.43 * HEATCAP_ELECTRON_SF * Tc * exp(-2.15 * (pow(Tc,2) - pow(temperature,2))/(Tc * temperature));
}

double calcElectronHC(double temperature){
    return temperature * HEATCAP_ELECTRON_SF;
}

double calcualteIc(double Ic0K, double Tc, double temperature){
    return Ic0K * pow(1.0-pow(temperature / Tc, 2),2);
}

double absFix( double input){
   if (input >= 0){
      return input;
   }
   else{
      return -1.0*input;
   }
}

// Simulation Function
bool isSC(double current, double temperature, double Ic0K , double Tc){
    return absFix(current) < calcualteIc(Ic0K, Tc, temperature) and temperature <= Tc;
}
void createHotspot(sSNSPD_X1 *opaque, double current){
    int PNR = 1;
    int hotspot_segments = (int)(opaque->sizehs * opaque->resolution / opaque->length);

    //("%i", hotspot_segments);
    double resistance = MINRES;

    for (int n = 0; n < PNR; ++n) {
        int start_index_hotspot = (1 + n)*opaque->resolution / (1 + PNR) - (hotspot_segments / 2);
        for (int i = start_index_hotspot; i < start_index_hotspot + hotspot_segments; ++i) {
            opaque->temperatures[i] = opaque->Ths;
            if (not isSC(current, opaque->temperatures[i], opaque->Ic0K, opaque->Tc)){
                resistance = opaque->Rsegment ++;
            }

        }
    }
    printf("%4.9f", resistance);
    opaque->resistance = resistance;

}

void diagonalSC(sSNSPD_X1 *opaque, int index, double dt, double* diagonal, double* off_diagonal, double* right_hand_side){
    double alpha = calcAlpha(opaque->temperatures[index]);
    double kapa = calcKapaSC(opaque->thickness, opaque->Rsheet, opaque->Tc, opaque->temperatures[index]);
    double heat_cap = calcElectronHCSC(opaque->Tc, opaque->temperatures[index]) + calcPhononHC(opaque->temperatures[index]);

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
    //cdef int i, it, il
    //cdef double m
    //cdef double resistance = 0
    double resistance = MINRES;
    double* diagonal = new double[opaque->resolution];
    double* off_diagonal = new double[opaque->resolution];
    double* right_hand_side = new double[opaque->resolution];

    diagonal[0] = diagonal[opaque->resolution - 1] = 1.0;
    right_hand_side[0] = right_hand_side[opaque->resolution - 1] = opaque->Tsub;
    off_diagonal[0] = off_diagonal[opaque->resolution - 1] = 0;

    //Get Diagonals for C-N matrix and total resistance

    for (int i = 1; i < opaque->resolution - 1; ++i) {
        if (isSC(current, opaque->temperatures[i], opaque->Ic0K, opaque->Tc)){
            diagonalSC(opaque, i, dt, diagonal, off_diagonal, right_hand_side);
        } else {
            diagonalNSC(opaque, i, dt, current, diagonal, off_diagonal, right_hand_side);
        }
    }
    double m;
    for (int it = 1; it < opaque->resolution; ++it) {
        m = off_diagonal[it] / diagonal[it-1];
        diagonal[it] = diagonal[it] - m * off_diagonal[it-1];
        right_hand_side[it] = right_hand_side[it] - m * right_hand_side[it-1];
    }


    opaque->temperatures[opaque->resolution - 1] = right_hand_side[opaque->resolution - 1] / diagonal[opaque->resolution -1];
    for (int il = opaque->resolution-2; il >=1; --il){
        opaque->temperatures[il] = (right_hand_side[il] - off_diagonal[il] *  opaque->temperatures[il + 1]) / diagonal[il];
        if (not isSC(current, opaque->temperatures[il], opaque->Ic0K, opaque->Tc)){
            resistance += opaque->Rsegment;
        }
    }

    opaque->resistance = resistance;
    free(diagonal);
    free(off_diagonal);
    free(right_hand_side);
}


//
sSNSPD_X1* initStates(union uData *data){
    int resolution = data[4].i;
    double Tsub = data[6].d;

    sSNSPD_X1* opaque = new sSNSPD_X1; // Use 'new' for C++ objects

    // If you still want to zero out the memory (though the constructor initializes)
    std::memset(opaque, 0, sizeof(sSNSPD_X1)); // C++ way to zero memory

    opaque->temperatures = new double[resolution]; // Use '->' to access member via pointer

    // Set all segments to
    for (int i = 0; i < resolution; ++i) {
        opaque->temperatures[i] = Tsub; // Example assignment
    }

    opaque->width      = data[1].d; // input parameter
    opaque->length     = data[2].d; // input parameter
    opaque->thickness  = data[3].d; // input parameter
    opaque->resolution = resolution; // input parameter
    opaque->Tc         = data[5].d; // input parameter
    opaque->Tsub       = Tsub; // input parameter
    opaque->Ic0K       = data[7].d; // input parameter
    opaque->Rsheet     = data[8].d; // input parameter
    opaque->hotspot    = data[9].b; // input parameter
    opaque->ths        = data[10].d; // input parameter
    opaque->Ths        = data[11].d; // input parameter
    opaque->sizehs     = data[12].d; // input parameter
    opaque->dlength    = opaque->length / (double)resolution;
    opaque->Rsegment   = opaque->Rsheet * opaque->dlength / opaque->width;
    return opaque;

}

extern "C" __declspec(dllexport) void snspd_x1(struct sSNSPD_X1 **opaque, double t, union uData *data)
{
   double current = data[ 0].d; // input
   double &ROUT = data[13].d; // output

   sSNSPD_X1 *inst = *opaque;

   if(!inst){
        inst = *opaque = initStates(data);
        ROUT=inst->resistance;
        inst->time = t;
        return;
        //printf("%f.8", inst->length);
        //inst->resolution;
   }
   double dt = t - inst->time;

   inst->time = t;
   if (inst->hotspot && t >= inst->ths){
        std::cout << "*click*";
        inst->hotspot = false;
        createHotspot(inst, current);
   }
   else{
        calcTotalResitance(inst, current, dt);
   }
   ROUT=inst->resistance;

}

extern "C" __declspec(dllexport) void Destroy(struct sSNSPD_X1 *inst)
{
    delete inst;
    inst = nullptr;
}
