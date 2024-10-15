#ifndef MATERIALS_H
#define MATERIALS_H

#include <QCoreApplication>
#include <QDebug>
#include <QMap>
#include <QPair>
#include <QString>
#include <QVector>
#include <dlfcn.h>

class Materials
{

public :

    // struct for complex refractive index of a material
    struct refractiveIndex {
        double realPart;
        double imaginaryPart;

        refractiveIndex(double re = 0, double im = 0) : realPart(re), imaginaryPart(im) {}
    };

    // refractive index function
    std::vector<refractiveIndex> RefractiveIndex(const std::vector<double>& energy, QString compound, double density)
    {
        std::vector<refractiveIndex> n;

        // Load the shared library
        void* handle = dlopen("/home/somayeh/xraylib-master/install/lib/libxrl.so.11", RTLD_LAZY);
        if (!handle) {
            qDebug() << "Cannot open the shared library ; "<< dlerror() << "\n";
            exit(EXIT_FAILURE);
        }

        // Load the function symbols
        typedef double (*Refractive_Index_ReFunc)(const char*, double, double, void*);
        typedef double (*Refractive_Index_ImFunc)(const char*, double, double, void*);

        Refractive_Index_ReFunc Refractive_Index_Re = (Refractive_Index_ReFunc) dlsym(handle, "Refractive_Index_Re");
        Refractive_Index_ImFunc Refractive_Index_Im = (Refractive_Index_ImFunc) dlsym(handle, "Refractive_Index_Im");

        const char* dlsym_error = dlerror();
        if (dlsym_error) {
            qDebug() << "Cannot load symbols: " << dlsym_error << "\n";
            dlclose(handle);
            exit(EXIT_FAILURE);
        }

        // Calculate the complex refractive index
        for (double E : energy) {
            refractiveIndex index;
            index.realPart = 1 - Refractive_Index_Re(compound.toStdString().c_str(), E, density, nullptr);
            index.imaginaryPart = Refractive_Index_Im(compound.toStdString().c_str(), E, density, nullptr);
            n.push_back(index);
        }

        // Close the library
        dlclose(handle);

        return n;
    }

    // Struct to hold the properties of materials
    struct MaterialProperties
    {
        double Z_A;      // Mean ratio of atomic number to mass
        double density;  // Density in g/cm^3
        QString formula; // Material formula
        QString name;    // Name of the material

        // MaterialProperties(double z_a = 0.0, double dens = 0.0, QString form = "", QString n = "")
        //     : Z_A(z_a), density(dens), formula(form), name(n) {}
    };

    // Initialize persistent data
    QVector<MaterialProperties> m_properties = {
        {0.99212, 8.375e-005, "H", "Hydrogen"},
        {0.49968, 1.663e-004, "He", "Helium"},
        {0.43221, 5.340e-001, "Li", "Lithium"},
        {0.44384, 1.848e+000, "Be", "Beryllium"},
        {0.46245, 2.370e+000, "B", "Boron"},
        {0.49954, 1.700e+000, "C", "Carbon"},
        {0.49976, 1.165e-003, "N", "Nitrogen"},
        {0.50002, 1.332e-003, "O", "Oxygen"},
        {0.47372, 1.580e-003, "F", "Fluorine"},
        {0.49555, 8.385e-004, "Ne", "Neon"},
        {0.47847, 9.710e-001, "Na", "Sodium"},
        {0.49373, 1.740e+000, "Mg", "Magnesium"},
        {0.48181, 2.699e+000, "Al", "Aluminum"},
        {0.49848, 2.330e+000, "Si", "Silicon"},
        {0.48428, 2.200e+000, "P", "Phosphorus"},
        {0.49897, 2.000e+000, "S", "Sulfur"},
        {0.47951, 2.995e-003, "Cl", "Chlorine"},
        {0.45059, 1.662e-003, "Ar", "Argon"},
        {0.48595, 8.620e-001, "K", "Potassium"},
        {0.49903, 1.550e+000, "Ca", "Calcium"},
        {0.46712, 2.989e+000, "Sc", "Scandium"},
        {0.45948, 4.540e+000, "Ti", "Titanium"},
        {0.45150, 6.110e+000, "V", "Vanadium"},
        {0.46157, 7.180e+000, "Cr", "Chromium"},
        {0.45506, 7.440e+000, "Mn", "Manganese"},
        {0.46556, 7.874e+000, "Fe", "Iron"},
        {0.45815, 8.900e+000, "Co", "Cobalt"},
        {0.47708, 8.902e+000, "Ni", "Nickel"},
        {0.45636, 8.960e+000, "Cu", "Copper"},
        {0.45879, 7.133e+000, "Zn", "Zinc"},
        {0.44462, 5.904e+000, "Ga", "Gallium"},
        {0.44071, 5.323e+000, "Ge", "Germanium"},
        {0.44046, 5.730e+000, "As", "Arsenic"},
        {0.43060, 4.500e+000, "Se", "Selenium"},
        {0.43803, 7.072e-003, "Br", "Bromine"},
        {0.42959, 3.478e-003, "Kr", "Krypton"},
        {0.43291, 1.532e+000, "Rb", "Rubidium"},
        {0.43369, 2.540e+000, "Sr", "Strontium"},
        {0.43867, 4.469e+000, "Y", "Yttrium"},
        {0.43848, 6.506e+000, "Zr", "Zirconium"},
        {0.44130, 8.570e+000, "Nb", "Niobium"},
        {0.43777, 1.022e+001, "Mo", "Molybdenum"},
        {0.43919, 1.150e+001, "Tc", "Technetium"},
        {0.43534, 1.241e+001, "Ru", "Ruthenium"},
        {0.43729, 1.241e+001, "Rh", "Rhodium"},
        {0.43225, 1.202e+001, "Pd", "Palladium"},
        {0.43572, 1.050e+001, "Ag", "Silver"},
        {0.42700, 8.650e+000, "Cd", "Cadmium"},
        {0.42676, 7.310e+000, "In", "Indium"},
        {0.42120, 7.310e+000, "Sn", "Tin"},
        {0.41889, 6.691e+000, "Sb", "Antimony"},
        {0.40752, 6.240e+000, "Te", "Tellurium"},
        {0.41764, 4.930e+000, "I", "Iodine"},
        {0.41130, 5.485e-003, "Xe", "Xenon"},
        {0.41383, 1.873e+000, "Cs", "Cesium"},
        {0.40779, 3.500e+000, "Ba", "Barium"},
        {0.41035, 6.154e+000, "La", "Lanthanum"},
        {0.41395, 6.657e+000, "Ce", "Cerium"},
        {0.41871, 6.710e+000, "Pr", "Praseodymium"},
        {0.41597, 6.900e+000, "Nd", "Neodymium"},
        {0.42094, 7.220e+000, "Pm", "Promethium"},
        {0.41234, 7.460e+000, "Sm", "Samarium"},
        {0.41457, 5.243e+000, "Eu", "Europium"},
        {0.40699, 7.900e+000, "Gd", "Gadolinium"},
        {0.40900, 8.229e+000, "Tb", "Terbium"},
        {0.40615, 8.550e+000, "Dy", "Dysprosium"},
        {0.40623, 8.795e+000, "Ho", "Holmium"},
        {0.40655, 9.066e+000, "Er", "Erbium"},
        {0.40844, 9.321e+000, "Tm", "Thulium"},
        {0.40453, 6.730e+000, "Yb", "Ytterbium"},
        {0.40579, 9.840e+000, "Lu", "Lutetium"},
        {0.40338, 1.331e+001, "Hf", "Hafnium"},
        {0.40343, 1.665e+001, "Ta", "Tantalum"},
        {0.40250, 1.930e+001, "W", "Tungsten"},
        {0.40278, 2.102e+001, "Re", "Rhenium"},
        {0.39958, 2.257e+001, "Os", "Osmium"},
        {0.40058, 2.242e+001, "Ir", "Iridium"},
        {0.39984, 2.145e+001, "Pt", "Platinum"},
        {0.40108, 1.932e+001, "Au", "Gold"},
        {0.39882, 1.355e+001, "Hg", "Mercury"},
        {0.39631, 1.172e+001, "Tl", "Thallium"},
        {0.39575, 1.135e+001, "Pb", "Lead"},
        {0.39717, 9.747e+000, "Bi", "Bismuth"},
        {0.40195, 9.320e+000, "Po", "Polonium"},
        {0.40479, 1.000e+001, "At", "Astatine"},
        {0.38736, 9.066e-003, "Rn", "Radon"},
        {0.39010, 1.000e+001, "Fr", "Francium"},
        {0.38934, 5.000e+000, "Ra", "Radium"},
        {0.39202, 1.007e+001, "Ac", "Actinium"},
        {0.38787, 1.172e+001, "Th", "Thorium"},
        {0.39388, 1.537e+001, "Pa", "Protactinium"},
        {0.38651, 1.895e+001, "U", "Uranium"},
        {0.39233, 2.030e+001, "Np", "Neptunium"},
        {0.38514, 1.984e+001, "Pu", "Plutonium"},
        {0.39085, 1.370e+001, "Am", "Americium"},
        {0.38855, 1.350e+001, "Cm", "Curium"},
        {0.39260, 1.400e+001, "Bk", "Berkelium"},
        {0.39031, 1.000e+001, "Cf", "Californium"},
        {0.39273, 8.840e+000, "Es", "Einsteinium"},
        {0.49897, 2.000e+000, "S", "Sulphur"},
        {0.41383, 1.873e+000, "Cs", "Caesium"},
        {0.49954, 1.700e+000, "C", "Graphite"},
        {0.49954, 3.513e+000, "C", "Diamond"},
        {0.49954, 2.000e+001, "C", "Carbon, Amorphous"},
        {0.55847, 9.200e-001, "H(0.119477)C(0.637240)N(0.007970)O(0.232333)Na(0.000500)Mg(0.000020)P(0.000160)S(0.000730)Cl(0.001190)K(0.000320)Ca(0.000020)Fe(0.000020)Zn(0.000020)", "Adipose tissue"},
        {0.49919, 1.205e-003, "C(0.000124)N(0.755267)O(0.231781)Ar(0.012827)", "Air"},
        {0.49038, 3.970e+000, "Al2O3", "Aluminum Oxide"},
        {0.55178, 1.100e+000, "C10H16O", "Amber"},
        {0.51129, 8.260e-004, "NH3", "Ammonia"},
        {0.53690, 1.024e+000, "C6H7N", "Aniline"},
        {0.52740, 1.283e+000, "C14H10", "Anthracene"},
        {0.52792, 1.250e+000, "H(0.057441)C(0.774591)O(0.167968)", "Bakelite"},
        {0.42207, 4.890e+000, "BaF2", "Barium Fluoride"},
        {0.44561, 4.500e+000, "BaSO4", "Barium Sulfate"},
        {0.47978, 3.010e+000, "BeO", "Beryllium Oxide"},
        {0.40961, 7.130e+000, "Bi12GeO20", "BGO"},
        {0.54995, 1.060e+000, "H(0.101866)C(0.100020)N(0.029640)O(0.759414)Na(0.001850)Mg(0.000040)Si(0.000030)P(0.000350)S(0.001850)Cl(0.002780)K(0.001630)Ca(0.000060)Fe(0.000460)Zn(0.000010)", "Blood"},
        {0.53010, 1.850e+000, "H(0.063984)C(0.278000)N(0.027000)O(0.410016)Mg(0.002000)P(0.070000)S(0.002000)Ca(0.147000)", "Compact Bone"},
        {0.52130, 1.850e+000, "H(0.047234)C(0.144330)N(0.041990)O(0.446096)Mg(0.002200)P(0.104970)S(0.003150)Ca(0.209930)Zn(0.000100)", "Cortical Bone"},
        {0.47058, 2.520e+000, "B4C", "Boron Carbide"},
        {0.48838, 1.812e+000, "B2O3", "Boron Oxide"},
        {0.55423, 1.030e+000, "H(0.110667)C(0.125420)N(0.013280)O(0.737723)Na(0.001840)Mg(0.000150)P(0.003540)S(0.001770)Cl(0.002360)K(0.003100)Ca(0.000090)Fe(0.000050)Zn(0.000010)", "Brain"},
        {0.55196, 1.020e+000, "H(0.106)C(0.332)N(0.03)O(0.527)Na(0.001)P(0.001)S(0.002)Cl(0.001)", "Breast"},
        {0.41665, 6.200e+000, "CdTe", "Cadmium Telluride"},
        {0.40749, 7.900e+000, "CdWO4", "Cadmium Tngstate"},
        {0.49955, 2.800e+000, "CaCO3", "Calcium Carbonate"},
        {0.48670, 3.180e+000, "CaF2", "Calcium Fluoride"},
        {0.49929, 3.300e+000, "O(0.285299)Ca(0.714701)", "Calcium Oxide"},
        {0.49950, 2.960e+000, "CaSO4", "Calcium Sulfate"},
        {0.40936, 6.062e+000, "CaWO4", "Calcium Tungstate"},
        {0.49989, 1.842e-003, "CO2", "Carbon Dioxide"},
        {0.42132, 4.115e+000, "CsF", "Cesium Fluoride"},
        {0.41569, 4.510e+000, "CsI", "Cesium Iodide"},
        {0.50274, 2.300e+000, "H(0.010000)C(0.001000)O(0.529107)Na(0.016000)Mg(0.002000)Al(0.033872)Si(0.337021)K(0.013000)Ca(0.044000)Fe(0.014000)", "Portland Concrete"},
        {0.50932, 2.300e+000, "H(0.0221)C(0.002484)O(0.57493)Na(0.015208)Mg(0.001266)Al(0.019953)Si(0.304627)K(0.010045)Ca(0.042951)Fe(0.006435)", "Concrete"},
        {0.45714, 3.350e+000,"H(0.003585)O(0.311622)Mg(0.001195)Al(0.004183)Si(0.010457)S(0.107858)Ca(0.050194)Fe(0.047505)Ba(0.4634)", "Barite Concrete"},
        {0.54877, 1.100e+000, "H(0.099269)C(0.193710)N(0.053270)O(0.653751)", "Eye lens"},
        {0.47592, 5.200e+000, "Fe2O3", "Iron(III) Oxide"},
        {0.46507, 7.150e+000, "FeB", "Ferroboride"},
        {0.47323, 5.700e+000, "FeO", "Iron(II) Oxide"},
        {0.53040, 1.420e+000, "C6H10O5", "Cellophane"},
        {0.53279, 1.200e+000, "C15H22O8", "CAB"},
        {0.51424, 1.490e+000, "C6H8N2O9", "Nitrocellulose"},
        {0.44247, 5.310e+000, "GaAs", "Gallium Arsenide"},
        {0.49707, 2.230e+000, "B(0.040064)O(0.539562)Na(0.028191)Al(0.011644)Si(0.377220)K(0.003321)", "Pyrex Glass"},
        {0.42101, 6.220e+000, "O(0.156453)Si(0.080866)Ti(0.008092)As(0.002651)Pb(0.751938)", "Lead Glass"},
        {0.49731, 2.400e+000, "O(0.459800)Na(0.096441)Si(0.336553)Ca(0.107205)", "Plate Glass"},
        {0.54292, 1.261e+000, "C3H8O3", "Glycerin"},
        {0.50039, 2.320e+000, "H4CaSO6", "Gypsum"},
        {0.40323, 9.530e+000, "PbO", "Lead Oxide"},
        {0.46262, 2.635e+000, "LiF", "Lithium Fluoride"},
        {0.41839, 3.494e+000, "LiI", "Lithium Iodide"},
        {0.46852, 2.013e+000, "Li2O", "Lithium Oxide"},
        {0.54965, 1.050e+000, "H(0.101278)C(0.102310)N(0.028650)O(0.757072)Na(0.001840)Mg(0.000730)P(0.000800)S(0.002250)Cl(0.002660)K(0.001940)Ca(0.000090)Fe(0.000370)Zn(0.000010)", "Lung"},
        {0.55512, 1.050e+000, "H(0.114318)C(0.655823)O(0.092183)Mg(0.134792)Ca(0.002883)", "M3 WAX"},
        {0.49814, 2.958e+000, "MgCO3", "Magnesium Carbonate"},
        {0.48153, 3.000e+000, "MgF2", "Magnesium Fluoride"},
        {0.49622, 3.580e+000, "MgO2", "Magnesium Oxide"},
        {0.49014, 2.530e+000, "MgB4O7", "Magnesium Tetraborate"},
        {0.52697, 1.200e+000, "C16H14O3", "Makrolon"},
        {0.40933, 6.360e+000, "HgI2", "Mercuric Iodide"},
        {0.45750, 8.60e+000, "CuZn5", "Brass"},
        {0.56479, 9.900e-001, "H(0.134040)C(0.777960)O(0.035020)Mg(0.038594)Ti(0.014386)", "MIX D WAX"},
        {0.53886, 1.000e+000, "H(0.081192)C(0.583442)N(0.017798)O(0.186381)Mg(0.130287)Cl(0.000900)", "MS20 Tissue Substitute"},
        {0.54938, 1.040e+000, "H(0.100637)C(0.107830)N(0.027680)O(0.754773)Na(0.000750)Mg(0.000190)P(0.001800)S(0.002410)Cl(0.000790)K(0.003020)Ca(0.000030)Fe(0.000040)Zn(0.000050)", "Skeletal muscle"},
        {0.55005, 1.040e+000, "H(0.101997)C(0.123000)N(0.035000)O(0.729003)Na(0.000800)Mg(0.000200)P(0.002000)S(0.005000)K(0.003000)", "Stirated muscle"},
        {0.54828, 1.110e+000, "H(0.098234)C(0.156214)N(0.035451)O(0.710100)", "Muscle-equivalent liquide, with sucrose"},
        {0.55014, 1.070e+000, "H(0.101969)C(0.120058)N(0.035451)O(0.742522)", "Muscle-equivalent liquide, without sucrose"},
        {0.53053, 1.145e+000, "C10H8", "Naphthalene"},
        {0.49985, 1.831e-003, "N2O", "Nitrous Oxide"},
        {0.55063, 1.080e+000, "H(0.103509)C(0.648415)N(0.099536)O(0.148539)", "Nylon, Du Pont Elvamide 8062"},
        {0.54790, 1.140e+000, "C12H22O2N2", "Nylon, Type 6 and Type 6/6"},
        {0.55236, 1.140e+000, "C8H15ON", "Nylon, Type 6/10"},
        {0.55649, 1.425e+000, "C11H21ON", "Nylon, Type 11 (Rilsan)"},
        {0.57778, 7.026e-001, "C8H18", "Liquid Octane"},
        {0.55149, 2.200e+000, "H(0.105)C(0.093)N(0.024)O(0.768)Na(0.002)P(0.002)S(0.002)Cl(0.002)K(0.002)", "Ovary"},
        {0.57275, 9.300e-001, "C25H52", "Paraffin"},
        {0.52767, 1.170e+000, "C3H3N", "Polyacrylonitrile"},
        {0.48558, 1.300e+000, "C17H18Cl2", "Polychlorostyrene"},
        {0.57034, 9.400e-001, "C2H4", "Polyethylene"},
        {0.52037, 1.400e+000, "C10H8O4", "Polyethylene Terephthalate(MYLAR)"},
        {0.53937, 1.190e+000, "C5H8O2", "Plexiglass"},
        {0.51183, 1.425e+000, "H2CO", "Polyoxymethylene"},
        {0.57034, 9.000e-001, "C3H6", "Polypropylene"},
        {0.53768, 1.060e+000, "CH", "Polysterene"},
        {0.47992, 2.200e+000, "C2F4", "Teflon"},
        {0.48081, 2.100e+000, "C2F3Cl", "PolytrifluorochloroEthylene"},
        {0.53432, 1.190e+000, "C4H6O2", "PVA"},
        {0.54480, 1.300e+000, "C2H4O", "PVOH"},
        {0.54537, 1.120e+000, "C8H13O2", "PVB"},
        {0.48710, 1.300e+000, "C2H3Cl", "PVC"},
        {0.49973, 1.760e+000, "C2H2F2", "Kynar"},
        {0.53984, 1.250e+000, "C6H9NO", "PVP"},
        {0.55785, 9.200e-001, "C5H8", "Latex"},
        {0.49930, 2.320e+000, "SiO2", "Silicon Oxide"},
        {0.42594, 6.010e+000, "AgI", "Silver Iodide"},
        {0.54932, 1.100e+000, "H(0.100588)C(0.228250)N(0.046420)O(0.619002)Na(0.000070)Mg(0.000060)P(0.000330)S(0.001590)Cl(0.002670)K(0.000850)Ca(0.000150)Fe(0.000010)Zn(0.000010)", "Skin"},
        {0.49062, 2.532e+000, "Na2CO3", "Sodium Carbonate"},
        {0.42697, 3.667e+000, "NaI", "Sodium Iodide"},
        {0.53170, 1.581e+000, "C12H22O11", "Sugar"},
        {0.55108, 1.040e+000, "H(0.104166)C(0.092270)N(0.019940)O(0.773884)Na(0.002260)Mg(0.000110)P(0.001250)S(0.001460)Cl(0.002440)K(0.002080)Ca(0.000100)Fe(0.000020)Zn(0.000020)", "Testes"},
        {0.55121, 1.000e+000, "H(0.104472)C(0.232190)N(0.024880)O(0.630238)Na(0.001130)Mg(0.000130)P(0.001330)S(0." "001990)Cl(0.001340)K(0.001990)Ca(0.000230)Fe(0.000050)Zn(0.000030)", "Soft Tissue"},
        {0.54975, 1.000e+000, "H(0.101172)C(0.111000)N(0.026000)O(0.761828)", "4-Component Soft Tissue"},
        {0.54993, 1.064e-003, "H(0.101869)C(0.456179)N(0.035172)O(0.406780)", "Tissue-equivalent gas (Methane-based)"},
        {0.55027, 1.826e-003, "H(0.102672)C(0.568940)N(0.035022)O(0.293366)", "Tissue-equivalent gas (Propane-based)"},
        {0.46528, 4.260e+000, "TiO2", "Titanium Oxide"},
        {0.40371, 2.400e+000, "WF6", "Tungstate HexaFluride"},
        {0.53284, 1.323e+000, "N2H4CO", "Urea"},
        {0.42266, 7.320e+000, "Gd2O2S", "CCD-Gadox"},
        {0.55509, 1.000e+000, "H2O", "Water"},
        {0.55509, 0.897e+000, "H2O", "Ice"},
        {0.55493, 0.93, "B(5)H(11.6)C(61.2)O(22.2)", "BPE"},
        {0.55012, 1.03, "B(10)H(11)C(58)O(21)", "BPE-10"},
        {0.46556, 7.86, "C(8.5E-4)P(1.5E-4)S(1.75E-4)Mn(0.003)Fe(0.996)", "Steel"}};

    // QMap<double, QString> atomicMass;

    // for (i = 0; i < properties.size(); i++) {
    //     atomicMass[properties[i].name] = properties[i].Z_A;
    // }

    // QPair<int, int> PhysProps() { return std::make_pair(properties.); }

    // Function to return physical properties of selected materials
    QVector<MaterialProperties> PhysProps(QVector<QString> names) {
        MaterialProperties result;
        QVector<MaterialProperties> result_vec;

        if (names.isEmpty()) {
            throw std::invalid_argument("No names provided");
        }

        for (const auto& Name : names) {
            bool found = false;
            for (int i = 0 ; i < m_properties.size() ; i++){
            // for (const auto& prop : properties) {
                if (Name == m_properties[i].formula || Name == m_properties[i].name) {
                    result = m_properties[i];
                    result_vec.push_back(result);
                    found = true;
                    break;
                }
            }
            if (!found) {
                throw std::invalid_argument("Material not found: " + Name.toStdString());
            }
        }

        return result_vec;
    }

};

#endif // MATERIALS_H
