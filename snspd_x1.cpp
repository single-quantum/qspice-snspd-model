//  cl /LD /std:c++17 /EHsc snspd_x1.cpp kernel32.lib

#include <math.h>
#include <execution>
#include <vector>
#include <numeric>
#include <algorithm> 

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
    std::vector<double> temperatures_curr;
    std::vector<double> temperatures_next;
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
    std::vector<double> diagonal;
    std::vector<double> off_diagonal;
    std::vector<double> right_hand_side;
    double  Tmax;
    double Lkin;
    double starttime;
    double Rmin;
    double TsubEr;
    double kapa_const;
    double kapa_sc_const;
    double r_const;
    double h_const;
    double g_sc_const;
    double g_const;
    double _1overTc;
    double _1overTcsquared;

    sSNSPD_X1(union uData *data) : resistance(data[15].d), time(0),
    width(data[1].d), length(data[2].d), thickness(data[3].d), resolution(data[4].i), Tsub(data[6].d),
    Tc(data[5].d),  Ic0K(data[7].d), Rsheet(data[8].d), hotspot(data[9].b), ths(data[10].d), 
    Ths(data[11].d), sizehs(data[12].d), photonnumber(data[13].i), Lkin(data[14].d), Rmin(data[15].d) {
        TsubEr = Tsub*(1+data[16].d);

        temperatures_curr.resize(resolution, Tsub);
        temperatures_next.resize(resolution, Tsub);
        diagonal.resize(resolution);
        off_diagonal.resize(resolution);
        right_hand_side.resize(resolution);

        const double Tc_cubed = Tc * Tc * Tc;

        GAMMA = data[19].d/Tc;
        A_const = 2.43*GAMMA*Tc;
        B_const = data[17].d/Tc_cubed;
        G_const = data[18].d/Tc_cubed;

        dlength    = length / (double)resolution;
        Rsegment   = Rsheet * dlength / width;
        Tmax       = Tsub;


        kapa_const = LORENTZ / (Rsheet * thickness);
        kapa_sc_const = kapa_const/ Tc;
        r_const = 1 / (2 * dlength* dlength);
        h_const = 1 / (2 * thickness);
        g_sc_const = Tsub / thickness;
        g_const = Rsheet/(width*width*thickness);

        _1overTc = 1/Tc;
        _1overTcsquared = _1overTc * _1overTc;
    }

};

// ... (DllMain, #undefs, and all calc* functions remain the same) ...
int __stdcall DllMain(void *module, unsigned int reason, void *reserved) { return 1; }

#undef CURRENT
#undef ROUT


double calcPow( double Tc,  double T, double _1overTc){
    return -2.15 * (Tc/T - T*_1overTc);
}

bool isSC(double current, double temperature, double Ic0K , double Tc, double _1overTcsquared){
    const double v = 1.0 - temperature * temperature * _1overTcsquared;
    return fabs(current) < Ic0K * v * v && temperature <= Tc;
}

void createHotspot(sSNSPD_X1 *opaque, double current){
    int PNR = opaque->photonnumber;
    int hotspot_segments = (int)(opaque->sizehs * opaque->resolution / opaque->length);
    double resistance = opaque->Rmin;
    std::vector<double>& T_curr = opaque->temperatures_curr;

    for (int n = 0; n < PNR; ++n) {
        int start_index_hotspot = (1 + n)*opaque->resolution / (1 + PNR) - (hotspot_segments / 2);
        for (int i = start_index_hotspot; i < start_index_hotspot + hotspot_segments; ++i) {
            T_curr[i] = opaque->Ths;
            if (!isSC(current, T_curr[i], opaque->Ic0K, opaque->Tc, opaque->_1overTcsquared)){
                resistance = resistance+ opaque->Rsegment;
            }
        }
    }
    opaque->resistance = resistance;
    opaque->Tmax = opaque->Ths;
}

void diagonalSC(sSNSPD_X1 *opaque, int index, double dt){
    const double& T_i =  opaque->temperatures_curr[index];
    const double& T_sum = opaque->temperatures_curr[index+1] + opaque->temperatures_curr[index-1];
    
    const double T_i_cubed = T_i * T_i * T_i;
    const double alpha = opaque->B_const * T_i_cubed;
    const double kapa = opaque->kapa_sc_const * T_i_cubed;
    const double heat_cap = opaque->A_const * exp(calcPow(opaque->Tc, T_i, opaque->_1overTc)) +  opaque->G_const * T_i_cubed;


    const double inv_HC = 1 / heat_cap;
    const double r = opaque->r_const*kapa * dt*inv_HC;
    const double h = opaque->h_const*alpha * dt*inv_HC;
    const double g = opaque->g_sc_const * alpha * dt*inv_HC;

    opaque->off_diagonal[index] = -r;
    opaque->diagonal[index]=  1 + h + 2 * r;
    opaque->right_hand_side[index] =  T_i  * (1 - h - 2 * r) + r * T_sum + g;
}

void  diagonalNSC(sSNSPD_X1 *opaque, int index, double dt, double current){
    const double& T_i =  opaque->temperatures_curr[index];
    const double& T_sum = opaque->temperatures_curr[index+1] + opaque->temperatures_curr[index-1];

    const double T_i_cubed = T_i * T_i * T_i;
    const double alpha = opaque->B_const * T_i_cubed;
    const double kapa = opaque->kapa_const * T_i;
    const double heat_cap = opaque->GAMMA * T_i + opaque->G_const * T_i_cubed;

    const double inv_HC = 1 / heat_cap;
    const double r = opaque->r_const*kapa * dt * inv_HC;
    const double h = opaque->h_const*alpha * dt * inv_HC;
    const double g = (opaque->g_sc_const * alpha  + current * current * opaque->g_const)*dt * inv_HC;

    opaque->off_diagonal[index] = -r;
    opaque->diagonal[index]=  1 + h + 2 * r;
    opaque->right_hand_side[index] = T_i  * (1 - h - 2 * r) + r * T_sum + g;
}

void calcTotalResitance(sSNSPD_X1 *opaque, double current, double dt){
    std::vector<double>& T_curr = opaque->temperatures_curr;
    std::vector<double>& T_next = opaque->temperatures_next;

    double resistance = opaque->Rmin;
    opaque->Tmax = opaque->Tsub;

    
    T_next[0] = T_next[opaque->resolution - 1] = opaque->Tsub;


    // Boundary conditions
    opaque->diagonal[0] = opaque->diagonal[opaque->resolution - 1] = 1.0;
    opaque->right_hand_side[0] = opaque->right_hand_side[opaque->resolution - 1] = opaque->Tsub;
    opaque->off_diagonal[0] = opaque->off_diagonal[opaque->resolution - 1] = 0;

    std::vector<int> indices(opaque->resolution - 2);
    std::iota(indices.begin(), indices.end(), 1);

    // Parallel construction of the matrix elements (safe)
    std::for_each(std::execution::par_unseq, indices.begin(), indices.end(),
                  [&](int i) {
                      if (isSC(current, T_curr[i], opaque->Ic0K, opaque->Tc, opaque->_1overTcsquared)){
                          diagonalSC(opaque, i, dt);
                      } else {
                          diagonalNSC(opaque, i, dt, current);
                      }
                  });

    // Forward Sweep (Sequential)
    for (int it = 1; it < opaque->resolution; ++it) {
        const double m = opaque->off_diagonal[it] / opaque->diagonal[it-1];
        opaque->diagonal[it] = opaque->diagonal[it] - m * opaque->off_diagonal[it-1];
        opaque->right_hand_side[it] = opaque->right_hand_side[it] - m * opaque->right_hand_side[it-1];
    }
    
    // Back Substitution (Sequential)
    T_next[opaque->resolution - 1] = opaque->right_hand_side[opaque->resolution - 1] / opaque->diagonal[opaque->resolution -1];
    
    for (int il = opaque->resolution-2; il >= 1; --il){
        T_next[il] = (opaque->right_hand_side[il] - opaque->off_diagonal[il] *  T_next[il + 1]) / opaque->diagonal[il];
        if (!isSC(current, T_next[il], opaque->Ic0K, opaque->Tc, opaque->_1overTcsquared)){
            resistance += opaque->Rsegment;
        }
        if (T_next[il] > opaque->Tmax){
            opaque->Tmax = T_next[il];
        }

    }

    T_curr.swap(T_next); 
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
   else if(inst->Tmax <= inst->TsubEr && isSC(CURRENT, inst->Tmax, inst->Ic0K, inst->Tc, inst->_1overTcsquared)){
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