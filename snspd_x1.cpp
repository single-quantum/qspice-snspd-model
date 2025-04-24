// Automatically generated C++ file on Thu Apr 24 13:03:17 2025
//
// To build with Digital Mars C++ Compiler:
//
//    dmc -mn -WD snspd_x1.cpp kernel32.lib

#include <malloc.h>
#include <math.h>
#include <iostream>
#include <cstring> // For bzero (though not the C++ way)


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
    double resitance;
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
    bool    hotspot;
    double  ths;
    double  Ths;
    double  sizehs;

    // Constructor (optional, but good practice in C++)
    sSNSPD_X1() : resitance(0.0), temperatures(nullptr) {}

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


// Simulation Function
void createHotspot(sSNSPD_X1 *opaque){
    int hotspot_segments = (int)(opaque->sizehs * opaque->resolution / opaque->length);
    float resistance = 0.0;

    for (int i = 0; i < opaque->resolution; ++i) {
        int start_index_hotspot = (1 + i)*opaque->resolution;
        for (int n = start_index_hotspot; n < start_index_hotspot + hotspot_segments; ++n) {
            int start_index_hotspot = (1 + i)*opaque->resolution;
            opaque->temperatures[i] = opaque->Ths;
            if (!){
                opaque->resistance += opaque->Rsheet;
            }

        }
    }

}



// 
sSNSPD_X1* initStates(union uData *data){
    int resolution = data[ 4].i;
    double Tsub = data[ 6].d;



    sSNSPD_X1* opaque = new sSNSPD_X1; // Use 'new' for C++ objects

    // If you still want to zero out the memory (though the constructor initializes)
    std::memset(opaque, 0, sizeof(sSNSPD_X1)); // C++ way to zero memory

    opaque->temperatures = new double[resolution]; // Use '->' to access member via pointer

    // Set all segments to 
    for (int i = 0; i < resolution; ++i) {
        opaque->temperatures[i] = Tsub; // Example assignment
    }

    std::cout << "First temperature: " << opaque->temperatures[0] << std::endl;

    opaque->width      = data[ 1].d; // input parameter
    opaque->length     = data[ 2].d; // input parameter
    opaque->thickness  = data[ 3].d; // input parameter
    opaque->resolution = resolution; // input parameter
    opaque->Tc         = data[ 5].d; // input parameter
    opaque->Tsub       = Tsub; // input parameter
    opaque->Ic0K       = data[ 7].d; // input parameter
    opaque->Rsheet     = data[ 8].d; // input parameter
    opaque->hotspot    = data[ 9].b; // input parameter
    opaque->ths        = data[10].d; // input parameter
    opaque->Ths        = data[11].d; // input parameter
    opaque->sizehs     = data[12].d; // input parameter

    return opaque;

}

extern "C" __declspec(dllexport) void snspd_x1(struct sSNSPD_X1 **opaque, double t, union uData *data)
{
   double  IN         = data[ 0].d; // input
   double &OUT        = data[13].d; // output

   if(!*opaque)
   {
        sSNSPD_X1* opaque = initStates(data);

   }
   struct sSNSPD_X1 *inst = *opaque;

// Implement module evaluation code here:

}

extern "C" __declspec(dllexport) void Destroy(struct sSNSPD_X1 *inst)
{
    delete inst;
    inst = nullptr;
   //free(inst);
}
