//    cl /LD /std:c++17 snspd_x1.cpp kernel32.lib

#include <malloc.h>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <execution>
#include <vector>
#include <numeric> // Added for std::iota
#include <algorithm> // Ensure std::fill is available
#include <cstring> // Added for memcpy

const double LORENTZ = 2.45e-8;

// ... (uData union remains the same) ...
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
    double B_const;
    double G_const;
    double A_const;
    double GAMMA;
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
    double Lkin;
    double starttime;
    double Rmin;
    double dt;
    double mod_counter;
    double dtmax;
    double TsubEr;

    // Constructor (optional, but good practice in C++)
    sSNSPD_X1(union uData *data) : resistance(data[15].d), temperatures(nullptr), time(0),
    width(data[1].d), length(data[2].d), thickness(data[3].d), resolution(data[4].i), Tsub(data[6].d),
    Tc(data[5].d),  Ic0K(data[7].d), Rsheet(data[8].d), hotspot(data[9].b), ths(data[10].d), 
    Ths(data[11].d), sizehs(data[12].d), photonnumber(data[13].i), Lkin(data[14].d), Rmin(data[15].d) {
        TsubEr = Tsub*(1+data[16].d);

        temperatures = new double[resolution]; // Use '->' to access member via pointer
        diagonal = new double[resolution];
        off_diagonal = new double[resolution];
        right_hand_side = new double[resolution];

        GAMMA = data[19].d/Tc;
        A_const = 2.43*GAMMA*Tc;
        B_const = data[17].d/pow(Tc,3);
        G_const = data[18].d/pow(Tc,3);

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

// ... (DllMain, #undefs, and all calc* functions remain the same) ...

int __stdcall DllMain(void *module, unsigned int reason, void *reserved) { return 1; }

#undef IN1
#undef IN2
#undef ROUT

double calcAlpha(double B,double temperature){
    return B * pow(temperature,3);
}

double calcKapa(double thickness, double Rsheet, double temperature){
    return ( LORENTZ * temperature ) / (Rsheet * thickness);
}

double calcKapaSC(double thickness, double Rsheet, double Tc, double temperature){
    return calcKapa(thickness, Rsheet, temperature) * temperature / Tc;
}

double calcPhononHC(double G, double temperature){
    return G * pow(temperature,3);
}

double calcPow( double Tc,  double T){
    return -2.15*Tc*(1-pow(T/Tc,2))/T;
}

double calcElectronHCSC(double Tc, double temperature, double A){
    return  A * exp(calcPow(Tc, temperature));
}

double calcElectronHC(double GAMMA, double temperature){
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

    for (int n = 0; n < PNR; ++n) {
        int start_index_hotspot = (1 + n)*opaque->resolution / (1 + PNR) - (hotspot_segments / 2);
        for (int i = start_index_hotspot; i < start_index_hotspot + hotspot_segments; ++i) {
            opaque->temperatures[i] = opaque->Ths;
            if (!isSC(current, opaque->temperatures[i], opaque->Ic0K, opaque->Tc)){
                resistance = resistance+ opaque->Rsegment;
            }
        }
    }
    opaque->resistance = resistance;
    opaque->Tmax = opaque->Ths;
}

void diagonalSC(sSNSPD_X1 *opaque, int index, double dt, double* diagonal, double* off_diagonal, double* right_hand_side){
    double alpha = calcAlpha(opaque->B_const, opaque->temperatures[index]);
    double kapa = calcKapaSC(opaque->thickness, opaque->Rsheet, opaque->Tc, opaque->temperatures[index]);
    double heat_cap = calcElectronHCSC(opaque->Tc, opaque->temperatures[index], opaque->A_const) + calcPhononHC(opaque->G_const, opaque->temperatures[index]);

    double r = kapa * dt / (2 * pow(opaque->dlength,2) * heat_cap);
    double h = alpha * dt / (2 * opaque->thickness * heat_cap);
    double g = dt / heat_cap * ( opaque->Tsub * alpha / opaque->thickness);

    off_diagonal[index] = -r;
    diagonal[index]=  1 + h + 2 * r;
    right_hand_side[index] = opaque->temperatures[index] * (1 - h - 2 * r) + r * (opaque->temperatures[index+1] + opaque->temperatures[index-1]) + g;
}

void  diagonalNSC(sSNSPD_X1 *opaque, int index, double dt, double current, double* diagonal, double* off_diagonal, double* right_hand_side){
    double alpha = calcAlpha(opaque->B_const, opaque->temperatures[index]);
    double kapa = calcKapa(opaque->thickness, opaque->Rsheet, opaque->temperatures[index]);
    double heat_cap = calcElectronHC(opaque->GAMMA, opaque->temperatures[index]) + calcPhononHC(opaque->G_const, opaque->temperatures[index]);

    double r = kapa * dt / (2 * pow(opaque->dlength,2) * heat_cap);
    double h = alpha * dt / (2 * opaque->thickness * heat_cap);
    double g = dt / heat_cap * ( opaque->Tsub * alpha  + pow(current / opaque->width, 2) * opaque->Rsheet) / opaque->thickness;

    off_diagonal[index] = -r;
    diagonal[index]=  1 + h + 2 * r;
    right_hand_side[index] = opaque->temperatures[index] * (1 - h - 2 * r) + r * (opaque->temperatures[index+1] + opaque->temperatures[index-1]) + g;
}

void calcTotalResitance(sSNSPD_X1 *opaque, double current, double dt){
    double resistance = opaque->Rmin;
    opaque->Tmax = opaque->Tsub;

    // Boundary conditions
    opaque->diagonal[0] = opaque->diagonal[opaque->resolution - 1] = 1.0;
    opaque->right_hand_side[0] = opaque->right_hand_side[opaque->resolution - 1] = opaque->Tsub;
    opaque->off_diagonal[0] = opaque->off_diagonal[opaque->resolution - 1] = 0;

    std::vector<int> indices(opaque->resolution - 2);
    std::iota(indices.begin(), indices.end(), 1);

    // Parallel construction of the matrix elements (safe)
    std::for_each(std::execution::par_unseq, indices.begin(), indices.end(),
                  [&](int i) {
                      if (isSC(current, opaque->temperatures[i], opaque->Ic0K, opaque->Tc)){
                          diagonalSC(opaque, i, dt, opaque->diagonal, opaque->off_diagonal, opaque->right_hand_side);
                      } else {
                          diagonalNSC(opaque, i, dt, current, opaque->diagonal, opaque->off_diagonal, opaque->right_hand_side);
                      }
                  });
    
    // The rest of the TDMA (Thomas Algorithm) must remain sequential!
    double m;
    double ai;
    double ci1;
    double* x;
    
    // Allocate the temporary array for the new temperatures
    x = new double[opaque->resolution];

    // Forward Sweep (Sequential)
    for (int it = 1; it < opaque->resolution; ++it) {
        ai = opaque->off_diagonal[it];
        ci1 = opaque->off_diagonal[it-1];
        m = ai / opaque->diagonal[it-1];
        opaque->diagonal[it] = opaque->diagonal[it] - m * ci1;
        opaque->right_hand_side[it] = opaque->right_hand_side[it] - m * opaque->right_hand_side[it-1];
    }
    
    // Back Substitution (Sequential)
    x[opaque->resolution - 1] = opaque->right_hand_side[opaque->resolution - 1] / opaque->diagonal[opaque->resolution -1];
    
    for (int il = opaque->resolution-2; il >= 1; --il){
        x[il] = (opaque->right_hand_side[il] - opaque->off_diagonal[il] *  x[il + 1]) / opaque->diagonal[il];
        if (!isSC(current, x[il], opaque->Ic0K, opaque->Tc)){
            resistance += opaque->Rsegment;
        }
        if (x[il] > opaque->Tmax){
            opaque->Tmax = x[il];
        }

    }

    // Set boundary condition for index 0
    x[0] = opaque->Tsub; 
    
    // Copy the new temperatures and clean up temporary array
    memcpy(opaque->temperatures, x, opaque->resolution * sizeof(double));
    delete[] x;

    opaque->resistance = resistance;
}


extern "C" __declspec(dllexport) void snspd_x1(struct sSNSPD_X1 **opaque, double t, union uData *data)
{
   double CURRENT = data[ 0].d; // input: current (I)
   double &SNSPRD_R = data[20].d; // output: resistance (R)
   sSNSPD_X1 *inst = *opaque;

   if(!inst){
        inst = *opaque = new sSNSPD_X1(data);
        inst->starttime = t;
   }

   if (t==inst->starttime) {
        inst->time = t;
        inst->resistance=inst->Rmin;
        SNSPRD_R = inst->resistance;
        return;
   }

   double dt = t - inst->time;
   inst->time = t;
   if (inst->hotspot && t >= inst->ths){
        inst->hotspot = false;
        createHotspot(inst, CURRENT);
   }
   else if(inst->Tmax <= inst->TsubEr && isSC(CURRENT, inst->Tmax, inst->Ic0K, inst->Tc)){
        inst->resistance = inst->Rmin;
        inst->Tmax = inst->Tsub;
   }
   else{
    calcTotalResitance(inst, CURRENT, dt);
   }
   SNSPRD_R = inst->resistance;
}

extern "C" __declspec(dllexport) void Destroy(struct sSNSPD_X1 *inst)
{
    delete inst;
    inst = nullptr;
}