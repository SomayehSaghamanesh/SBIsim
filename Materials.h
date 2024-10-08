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
        {0.48181, 2.699e+000, "Al", "Aluminium"},
        {0.49897, 2.000e+000, "S", "Sulphur"},
        {0.41383, 1.873e+000, "Cs", "Caesium"},
        {0.49954, 1.700e+000, "C", "Graphite"},
        {0.49954, 3.513e+000, "C", "Diamond"},
        {0.49954, 2.000e+001, "C", "Carbon, Amorphous"},
        {0.54903, 1.127e+000, "H(0.101327)C(0.775501)N(0.035057)O(0.052316)F(0.017422)Ca(0.018378)", "A-150"},
        {0.55097, 7.899e-001, "C3H6O", "ACETONE"},
        {0.53768, 1.097e-003, "C2H2", "ACETYLENE"},
        {0.51803, 1.350e+000, "C5H5N5", "ADENINE"},
        {0.55847, 9.200e-001, "H(0.119477)C(0.637240)N(0.007970)O(0.232333)Na(0.000500)Mg(0.000020)P(0.000160)S(0.000730)Cl(0.001190)K(0.000320)Ca(0.000020)Fe(0.000020)Zn(0.000020)", "ADIPOSE TISSUE"},
        {0.49919, 1.205e-003, "C(0.000124)N(0.755267)O(0.231781)Ar(0.012827)", "Air"},
        {0.53876, 1.420e+000, "C3H7NO2", "ALANINE"},
        {0.49038, 3.970e+000, "Al2O3", "ALUMINUM OXIDE"},
        {0.55178, 1.100e+000, "C10H16O", "AMBER"},
        {0.51129, 8.260e-004, "NH3", "AMMONIA"},
        {0.53690, 1.024e+000, "C6H7N", "ANILINE"},
        {0.52740, 1.283e+000, "C14H10", "ANTHRACENE"},
        {0.52740, 1.450e+000, "H(0.065471)C(0.536945)N(0.021500)O(0.032085)F(0.167411)Ca(0.176589)", "B-100"},
        {0.52792, 1.250e+000, "H(0.057441)C(0.774591)O(0.167968)", "BAKELITE"},
        {0.42207, 4.890e+000, "BaF2", "BARIUM FLUORIDE"},
        {0.44561, 4.500e+000, "BaSO4", "BARIUM SULFATE"},
        {0.53768, 8.787e-001, "C6H6", "BENZENE"},
        {0.47978, 3.010e+000, "BeO", "BERYLLIUM OXIDE"},
        {0.40961, 7.130e+000, "Bi12GeO20", "BGO"},
        {0.54995, 1.060e+000, "H(0.101866)C(0.100020)N(0.029640)O(0.759414)Na(0.001850)Mg(0.000040)Si(0.000030)P(0.000350)S(0.001850)Cl(0.002780)K(0.001630)Ca(0.000060)Fe(0.000460)Zn(0.000010)", "Blood"},
        {0.53010, 1.850e+000, "H(0.063984)C(0.278000)N(0.027000)O(0.410016)Mg(0.002000)P(0.070000)S(0.002000)Ca(0.147000)", "COMPACT BONE"},
        {0.52130, 1.850e+000, "H(0.047234)C(0.144330)N(0.041990)O(0.446096)Mg(0.002200)P(0.104970)S(0.003150)Ca(0.209930)Zn(0.000100)", "CORTICAL BONE"},
        {0.47058, 2.520e+000, "B4C", "BORON CARBIDE"},
        {0.48838, 1.812e+000, "B2O3", "BORON OXIDE"},
        {0.55423, 1.030e+000, "H(0.110667)C(0.125420)N(0.013280)O(0.737723)Na(0.001840)Mg(0.000150)P(0.003540)S(0.001770)Cl(0.002360)K(0.003100)Ca(0.000090)Fe(0.000050)Zn(0.000010)", "Brain"},
        {0.55196, 1.020e+000, "H(0.106)C(0.332)N(0.03)O(0.527)Na(0.001)P(0.001)S(0.002)Cl(0.001)", "Breast"},
        {0.58497, 2.493e-003, "C4H10", "BUTANE"},
        {0.56663, 8.098e-001, "C4H10O", "n-BUTANOL"},
        {0.49969, 1.760e+000, "H(0.024680)C(0.501610)O(0.004527)F(0.465209)Si(0.003973)", "C-552"},
        {0.41665, 6.200e+000, "CdTe", "CADMIUM TELLURIDE"},
        {0.40749, 7.900e+000, "CdWO4", "CWO"},
        {0.49955, 2.800e+000, "CaCO3", "CALCIUM CARBONATE"},
        {0.48670, 3.180e+000, "CaF2", "CALCIUM FLUORIDE"},
        {0.49929, 3.300e+000, "O(0.285299)Ca(0.714701)", "CALCIUM OXIDE"},
        {0.49950, 2.960e+000, "CaSO4", "CALCIUM SULFATE"},
        {0.40936, 6.062e+000, "CaWO4", "CALCIUM TUNGSTATE"},
        {0.49989, 1.842e-003, "CO2", "CARBON DIOXIDE"},
        {0.48107, 1.594e+000, "CCl4", "CARBON TETRACHLORIDE"},
        {0.53040, 1.420e+000, "C6H10O5", "CELLOPHANE"},
        {0.53279, 1.200e+000, "C15H22O8", "CAB"},
        {0.51424, 1.490e+000, "C6H8N2O9", "Nitrocellulose"},
        {0.55278, 1.030e+000, "H(0.107596)N(0.000800)O(0.874976)S(0.014627)Ce(0.002001)", "CERIC SULFATE DOSIMETER SOLUTION"},
        {0.42132, 4.115e+000, "CsF", "CESIUM FLUORIDE"},
        {0.41569, 4.510e+000, "CsI", "CESIUM IODIDE"},
        {0.48734, 1.106e+000, "C6H5Cl", "CHLOROBENZENE"},
        {0.51498, 1.483e+000, "CHCl3", "CHLOROFORM"},
        {0.50274, 2.300e+000, "H(0.010000)C(0.001000)O(0.529107)Na(0.016000)Mg(0.002000)Al(0.033872)Si(0.337021)K(0.013000)Ca(0.044000)Fe(0.014000)", "Portland Concrete"},
        {0.50932, 2.300e+000, "H(0.0221)C(0.002484)O(0.57493)Na(0.015208)Mg(0.001266)Al(0.019953)Si(0.304627)K(0.010045)Ca(0.042951)Fe(0.006435)", "Concrete"},
        {0.45714, 3.350e+000,"H(0.003585)O(0.311622)Mg(0.001195)Al(0.004183)Si(0.010457)S(0.107858)Ca(0.050194)Fe(0.047505)Ba(0.4634)", "Barite Concrete"},
        {0.57034, 7.790e-001, "C6H12", "CYCLOHEXANE"},
        {0.49098, 1.305e+000, "C6H4Cl2", "ortho-dichlorobenzene"},
        {0.48616, 1.220e+000, "C4H8Cl2O", "DICHLORODIETHYL ETHER"},
        {0.48853, 1.235e+000, "C2H4Cl2", "1, 2-DICHLOROETHANE"},
        {0.56663, 7.138e-001, "C4H10O", "DIETHYL ETHER"},
        {0.54724, 9.487e-001, "C3H7NO", "N,N-DIMETHYL FORMAMIDE"},
        {0.53757, 1.101e+000, "C2H6OS", "DIMETHYL SULFOXIDE"},
        {0.59862, 1.253e-003, "C2H6", "ETHANE"},
        {0.56437, 7.893e-001, "C2H6O", "Ethanol"},
        {0.54405, 1.130e+000, "C12H22O5", "ETHYL CELLULOSE"},
        {0.57034, 1.175e-003, "", "ETHYLENE"},
        {0.54877, 1.100e+000, "H(0.099269)C(0.193710)N(0.053270)O(0.653751)", "EYE LENS"},
        {0.47592, 5.200e+000, "Fe2O3", "Iron(III) oxide"},
        {0.46507, 7.150e+000, "FeB", "FERROBORIDE"},
        {0.47323, 5.700e+000, "FeO", "Iron(II) oxide"},
        {0.58497, 2.493e-003, "C4H10", "BUTANE"},
        {0.56663, 8.098e-001, "C4H10O", "n-BUTANOL"},
        {0.49969, 1.760e+000, "H(0.024680)C(0.501610)O(0.004527)F(0.465209)Si(0.003973)", "C-552"},
        {0.41665, 6.200e+000, "CdTe", "CADMIUM TELLURIDE"},
        {0.40749, 7.900e+000, "CdWO4", "CWO"},
        {0.49955, 2.800e+000, "CaCO3", "CALCIUM CARBONATE"},
        {0.48670, 3.180e+000, "CaF2", "CALCIUM FLUORIDE"},
        {0.49929, 3.300e+000, "O(0.285299)Ca(0.714701)", "CALCIUM OXIDE"},
        {0.49950, 2.960e+000, "CaSO4", "CALCIUM SULFATE"},
        {0.40936, 6.062e+000, "CaWO4", "CALCIUM TUNGSTATE"},
        {0.49989, 1.842e-003, "CO2", "CARBON DIOXIDE"},
        {0.48107, 1.594e+000, "CCl4", "CARBON TETRACHLORIDE"},
        {0.53040, 1.420e+000, "C6H10O5", "CELLOPHANE"},
        {0.53279, 1.200e+000, "C15H22O8", "CAB"},
        {0.51424, 1.490e+000, "C6H8N2O9", "Nitrocellulose"},
        {0.55278, 1.030e+000, "H(0.107596)N(0.000800)O(0.874976)S(0.014627)Ce(0.002001)", "CERIC SULFATE DOSIMETER SOLUTION"},
        {0.42132, 4.115e+000, "CsF", "CESIUM FLUORIDE"},
        {0.41569, 4.510e+000, "CsI", "CESIUM IODIDE"},
        {0.48734, 1.106e+000, "C6H5Cl", "CHLOROBENZENE"},
        {0.51498, 1.483e+000, "CHCl3", "CHLOROFORM"},
        {0.50274, 2.300e+000, "H(0.010000)C(0.001000)O(0.529107)Na(0.016000)Mg(0.002000)Al(0.033872)Si(0.337021)K(0.013000)Ca(0.044000)Fe(0.014000)", "Portland Concrete"},
        {0.50932, 2.300e+000, "H(0.0221)C(0.002484)O(0.57493)Na(0.015208)Mg(0.001266)Al(0.019953)Si(0.304627)K(0.010045)Ca(0.042951)Fe(0.006435)", "Concrete"},
        {0.45714, 3.350e+000, "H(0.003585)O(0.311622)Mg(0.001195)Al(0.004183)Si(0.010457)S(0.107858)Ca(0.050194)Fe(0.047505)Ba(0.4634)", "Barite Concrete"},
        {0.57034, 7.790e-001, "C6H12", "CYCLOHEXANE"},
        {0.49098, 1.305e+000, "C6H4Cl2", "ortho-dichlorobenzene"},
        {0.48616, 1.220e+000, "C4H8Cl2O", "DICHLORODIETHYL ETHER"},
        {0.48853, 1.235e+000, "C2H4Cl2", "1, 2-DICHLOROETHANE"},
        {0.56663, 7.138e-001, "C4H10O", "DIETHYL ETHER"},
        {0.54724, 9.487e-001, "C3H7NO", "N,N-DIMETHYL FORMAMIDE"},
        {0.55328, 1.024e+000, "H(0.108259)N(0.000027)O(0.878636)Na(0.000022)S(0.012968)Cl(0.000034)Fe(0.000054)", "FERROUS SULFATE DOSIMETER SOLUTION"},
        {0.47968, 1.120e+000, "CCl2F2", "FREON-12"},
        {0.44801, 1.800e+000, "CBr2F2", "FREON-12B2"},
        {0.47866, 9.500e-001, "CClF3", "FREON-13"},
        {0.45665, 1.500e+000, "CBrF3", "FREON-13B1"},
        {0.42262, 1.800e+000, "CIF3", "FREON-13I1"},
        {0.42266, 7.440e+000, "Gd2O2S", "GADOLINIUM OXYSULFIDE"},
        {0.44247, 5.310e+000, "GaAs", "GALLIUM ARSENIDE"},
        {0.53973, 1.291e+000, "H(0.081180)C(0.416060)N(0.111240)O(0.380640)S(0.010880)", "GEL IN PHOTOGRAPHIC EMULSION"},
        {0.49707, 2.230e+000, "B(0.040064)O(0.539562)Na(0.028191)Al(0.011644)Si(0.377220)K(0.003321)", "Pyrex Glass"},
        {0.42101, 6.220e+000, "O(0.156453)Si(0.080866)Ti(0.008092)As(0.002651)Pb(0.751938)", "Lead Glass"},
        {0.49731, 2.400e+000, "O(0.459800)Na(0.096441)Si(0.336553)Ca(0.107205)", "Plate Glass"},
        {0.53489, 1.540e+000, "C6H14O7", "GLUCOSE"},
        {0.53371, 1.460e+000, "C5H10N2O3", "GLUTAMINE"},
        {0.54292, 1.261e+000, "C3H8O3", "Glycerin"},
        {0.51612, 1.580e+000, "C5H5N5O", "GUANINE"},
        {0.50039, 2.320e+000, "H4CaSO6", "GYPSUM"},
        {0.57882, 6.838e-001, "C7H16", "HEPTANE"},
        {0.58020, 6.603e-001, "C6H14", "HEXANE"},
        {0.51264, 1.420e+000, "C22H10N2O5", "KAPTON POLYIMIDE FILM"},
        {0.42588, 6.280e+000, "LaOBr", "LANTHANUM OXYBROMIDE"},
        {0.42706, 5.860e+000, "La2OS", "LANTHANUM OXYSULFIDE"},
        {0.40323, 9.530e+000, "PbO", "LEAD OXIDE"},
        {0.50052, 1.178e+000, "LiNH2", "LITHIUM AMIDE"},
        {0.48720, 2.110e+000, "Li2CO3", "LITHIUM CARBONATE"},
        {0.46262, 2.635e+000, "LiF", "LITHIUM FLUORIDE"},
        {0.50321, 8.200e-001, "LiH", "LITHIUM HYDRIDE"},
        {0.41839, 3.494e+000, "LiI", "LITHIUM IODIDE"},
        {0.46852, 2.013e+000, "Li2O", "LITHIUM OXIDE"},
        {0.48487, 2.440e+000, "Li2B4O7", "LITHIUM TETRABORATE"},
        {0.54965, 1.050e+000, "H(0.101278)C(0.102310)N(0.028650)O(0.757072)Na(0.001840)Mg(0.000730)P(0.000800)S(0.002250)Cl(0.002660)K(0.001940)Ca(0.000090)Fe(0.000370)Zn(0.000010)", "Lung"},
        {0.55512, 1.050e+000, "H(0.114318)C(0.655823)O(0.092183)Mg(0.134792)Ca(0.002883)", "M3 WAX"},
        {0.49814, 2.958e+000, "MgCO3", "MAGNESIUM CARBONATE"},
        {0.48153, 3.000e+000, "MgF2", "MAGNESIUM FLUORIDE"},
        {0.49622, 3.580e+000, "MgF2", "MAGNESIUM OXIDE"},
        {0.49014, 2.530e+000, "MgB4O7", "MAGNESIUM TETRABORATE"},
        {0.40933, 6.360e+000, "HgI2", "MERCURIC IODIDE"},
        {0.45750, 8.60e+000, "CuZn5", "messing"},
        {0.62334, 6.672e-004, "CH4", "METHANE"},
        {0.56176, 7.914e-001, "CH4O", "METHANOL"},
        {0.56479, 9.900e-001, "H(0.134040)C(0.777960)O(0.035020)Mg(0.038594)Ti(0.014386)", "MIX D WAX"},
        {0.53886, 1.000e+000, "H(0.081192)C(0.583442)N(0.017798)O(0.186381)Mg(0.130287)Cl(0.000900)", "MS20 TISSUE SUBSTITUTE"},
        {0.54938, 1.040e+000, "H(0.100637)C(0.107830)N(0.027680)O(0.754773)Na(0.000750)Mg(0.000190)P(0.001800)S(0.002410)Cl(0.000790)K(0.003020)Ca(0.000030)Fe(0.000040)Zn(0.000050)", "SKELETAL MUSCLE"},
        {0.55005, 1.040e+000, "H(0.101997)C(0.123000)N(0.035000)O(0.729003)Na(0.000800)Mg(0.000200)P(0.002000)S(0.005000)K(0.003000)", "STRIATED MUSCLE"},
        {0.54828, 1.110e+000, "H(0.098234)C(0.156214)N(0.035451)O(0.710100)", "MUSCLE-EQUIVALENT LIQUID, WITH SUCROSE"},
        {0.55014, 1.070e+000, "H(0.101969)C(0.120058)N(0.035451)O(0.742522)", "MUSCLE-EQUIVALENT LIQUID, WITHOUT SUCROSE"},
        {0.53053, 1.145e+000, "C10H8", "NAPHTHALENE"},
        {0.51986, 1.199e+000, "C6H5NO2", "NITROBENZENE"},
        {0.49985, 1.831e-003, "N2O", "NITROUS OXIDE"},
        {0.55063, 1.080e+000, "H(0.103509)C(0.648415)N(0.099536)O(0.148539)", "NYLON, DU PONT ELVAMIDE 8062"},
        {0.54790, 1.140e+000, "C12H22O2N2", "NYLON, TYPE 6 AND TYPE 6/6"},
        {0.55236, 1.140e+000, "C8H15ON", "NYLON, TYPE 6/10"},
        {0.55649, 1.425e+000, "C11H21ON", "NYLON, TYPE 11 (RILSAN)"},
        {0.57778, 7.026e-001, "C8H18", "Liquid Octane"},
        {0.55149, 2.200e+000, "H(0.105)C(0.093)N(0.024)O(0.768)Na(0.002)P(0.002)S(0.002)Cl(0.002)K(0.002)", "Ovary"},
        {0.57275, 9.300e-001, "C25H52", "PARAFFIN"},
        {0.58212, 6.262e-001, "C5H12", "PENTANE"},
        {0.48176, 3.820e+000, "H(0.0305)C(0.2107)N(0.0721)O(0.1632)Br(0.2228)Ag(0.3007)", "Kodak Photo Emulsion"},
        {0.45453, 1.030e+000, "H(0.0141)C(0.072261)N(0.01932)O(0.066101)S(0.00189)Br(0.349104)Ag(0.474105)I(0.00312)", "Nuclear Photo Emulsion"},
        {0.45453, 3.815e+000, "H(0.0141)C(0.072261)N(0.01932)O(0.066101)S(0.00189)Br(0.349103)Ag(0.474105)I(0.00312)", "PHOTOGRAPHIC EMULSION"},
        {0.54141, 1.032e+000, "C9H10", "PLASTIC SCINTILLATOR (VINYLTOLUENE BASED)"},
        {0.38879, 1.146e+001, "PuO2", "PLUTONIUM DIOXIDE"},
        {0.52767, 1.170e+000, "C3H3N", "POLYACRYLONITRILE"},
        {0.52697, 1.200e+000, "C16H14O3", "MAKROLON"},
        {0.48558, 1.300e+000, "C17H18Cl2", "POLYCHLOROSTYRENE"},
        {0.57034, 9.400e-001, "C2H4", "POLYETHYLENE"},
        {0.52037, 1.400e+000, "C10H8O4", "POLYETHYLENE TEREPHTHALATE (MYLAR)"},
        {0.53937, 1.190e+000, "C5H8O2", "PLEXIGLASS"},
        {0.51183, 1.425e+000, "H2CO", "POLYOXYMETHYLENE"},
        {0.57034, 9.000e-001, "C3H6", "POLYPROPYLENE"},
        {0.53768, 1.060e+000, "CH", "POLYSTYRENE"},
        {0.47992, 2.200e+000, "C2F4", "Teflon"},
        {0.48081, 2.100e+000, "C2F3Cl", "POLYTRIFLUOROCHLOROETHYLENE"},
        {0.53432, 1.190e+000, "C4H6O2", "PVA"},
        {0.54480, 1.300e+000, "C2H4O", "PVOH"},
        {0.54537, 1.120e+000, "C8H13O2", "PVB"},
        {0.48710, 1.300e+000, "C2H3Cl", "PVC"},
        {0.49513, 1.700e+000, "CHCl", "SARAN"},
        {0.49973, 1.760e+000, "C2H2F2", "Kynar"},
        {0.53984, 1.250e+000, "C6H9NO", "PVP"},
        {0.43373, 3.130e+000, "KI", "POTASSIUM IODIDE"},
        {0.48834, 2.320e+000, "K2O", "POTASSIUM OXIDE"},
        {0.58962, 1.879e-003, "C3H8", "PROPANE"},
        {0.58962, 4.300e-001, "C3H8", "LIQUID PROPANE"},
        {0.56576, 8.035e-001, "C3H8O", "1-Propanol"},
        {0.53097, 9.819e-001, "C5H5N", "PYRIDINE"},
        {0.57034, 9.200e-001, "C4H8", "BUTYL"},
        {0.55785, 9.200e-001, "C5H8", "Latex"},
        {0.48605, 1.230e+000, "C4H5Cl", "NEOPRENE"},
        {0.49930, 2.320e+000, "SiO2", "SILICON DIOXIDE"},
        {0.43670, 6.473e+000, "AgBr", "SILVER BROMIDE"},
        {0.44655, 5.560e+000, "AgCl", "SILVER CHLORIDE"},
        {0.43663, 6.470e+000, "Br(0.422895)Ag(0.573748)I(0.003357)", "SILVER HALIDES IN PHOTOGRAPHIC EMULSION"},
        {0.42594, 6.010e+000, "AgI", "SILVER IODIDE"},
        {0.54932, 1.100e+000, "H(0.100588)C(0.228250)N(0.046420)O(0.619002)Na(0.000070)Mg(0.000060)P(0.000330)S(0.001590)Cl(0.002670)K(0.000850)Ca(0.000150)Fe(0.000010)Zn(0.000010)", "Skin"},
        {0.49062, 2.532e+000, "Na2CO3", "SODIUM CARBONATE"},
        {0.42697, 3.667e+000, "NaI", "SODIUM IODIDE"},
        {0.48403, 2.270e+000, "Na2O", "SODIUM MONOXIDE"},
        {0.49415, 2.261e+000, "NaNO3", "SODIUM NITRATE"},
        {0.53260, 9.707e-001, "C14H12", "STILBENE"},
        {0.53170, 1.581e+000, "C12H22O11", "Sugar"},
        {0.52148, 1.234e+000, "C18H10", "TERPHENYL"},
        {0.55108, 1.040e+000, "H(0.104166)C(0.092270)N(0.019940)O(0.773884)Na(0.002260)Mg(0.000110)P(0.001250)S(0.001460)Cl(0.002440)K(0.002080)Ca(0.000100)Fe(0.000020)Zn(0.000020)", "TESTES"},
        {0.48241, 1.625e+000, "C2Cl4", "PCE"},
        {0.40861, 7.004e+000, "TlCl", "THALLIUM CHLORIDE"},
        {0.55121, 1.000e+000, "H(0.104472)C(0.232190)N(0.024880)O(0.630238)Na(0.001130)Mg(0.000130)P(0.001330)S(0." "001990)Cl(0.001340)K(0.001990)Ca(0.000230)Fe(0.000050)Zn(0.000030)", "Soft Tissue"},
        {0.54975, 1.000e+000, "H(0.101172)C(0.111000)N(0.026000)O(0.761828)", "4 component soft tissue"},
        {0.54993, 1.064e-003, "H(0.101869)C(0.456179)N(0.035172)O(0.406780)", "TISSUE-EQUIVALENT GAS (METHANE BASED)"},
        {0.55027, 1.826e-003, "H(0.102672)C(0.568940)N(0.035022)O(0.293366)", "TISSUE-EQUIVALENT GAS (PROPANE BASED)"},
        {0.46528, 4.260e+000, "TiO2", "TITANIUM DIOXIDE"},
        {0.54265, 8.669e-001, "C7H8", "TOLUENE"},
        {0.48710, 1.460e+000, "Cl3C2H", "TRICHLOROETHYLENE"},
        {0.52404, 1.070e+000, "C6H15PO4", "TRIETHYL PHOSPHATE"},
        {0.40371, 2.400e+000, "WF6", "TUNGSTEN HEXAFLUORIDE"},
        {0.39687, 1.128e+001, "UC2", "URANIUM DICARBIDE"},
        {0.39194, 1.363e+001, "UC", "URANIUM MONOCARBIDE"},
        {0.39996, 1.096e+001, "UO2", "URANIUM OXIDE"},
        {0.53284, 1.323e+000, "N2H4CO", "UREA"},
        {0.54632, 1.230e+000, "C5H11NO2", "VALINE"},
        {0.52772, 1.800e+000, "C5H2F8", "Viton"},
        {0.42266, 7.320e+000, "Gd2O2S", "CCD-Gadox"},
        {0.41569, 4.510e+000, "CsI", "Imager-CsI"},
        {0.55509, 1.000e+000, "H2O", "Water"},
        {0.55509, 7.562e-004, "H2O", "Steam"},
        {0.55509, 0.897e+000, "H2O", "Ice"},
        {0.54631, 8.700e-001, "C8H10", "XYLENE"},
        {0.38930, 2.055, "U3O8", "Yellowcake"},
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
