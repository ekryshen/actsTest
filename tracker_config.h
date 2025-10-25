#ifndef tracker_config
#define tracker_config

double bz = 0.5;    // T

//std::vector<double> positions{{210, 232.5, 255, 277.5, 300}}; // z-positions of sensitive layers in cm
//std::vector<double> positions{{210, 211, 232.5, 233.5, 255, 256, 277.5, 278.5, 300}}; // double layers for strip-like simulations
std::vector<double> positions{{210, 210.8, 211.6, 212., 212.4, 213.2, 214.,
                               232, 232.8, 233.6, 234., 234.4, 235.2, 236.,
                               254, 254.8, 255.6, 256., 256.4, 257.2, 258.,
                               276, 276.8, 277.6, 278., 278.4, 279.2, 280.,
                               298, 298.8, 299.6, 300., 300.4, 301.2, 302.
}}; // cm

std::vector<double> layerAngle {{0., 0., 0., 0., 0.0041, 0.0041, 0.0041,
                                 0., 0., 0., 0., 0.0041, 0.0041, 0.0041,
                                 0., 0., 0., 0., 0.0041, 0.0041, 0.0041,
                                 0., 0., 0., 0., 0.0041, 0.0041, 0.0041,
                                 0., 0., 0., 0., 0.0041, 0.0041, 0.0041
}}; // rad

std::vector<double> layerRMin {{61., 61., 61., 61., 61., 61., 61.,
                                67., 67., 67., 67., 67., 67., 67.,
                                73., 73., 73., 73., 73., 73., 73.,
                                80., 80., 80., 80., 80., 80., 80.,
                                86., 86., 86., 86., 86., 86., 86.
}}; // cm

std::vector<double> layerRMax{{  93., 93., 93., 93., 93., 93., 93.,
                                103.,103.,103.,103.,103.,103.,103.,
                                113.,113.,113.,113.,113.,113.,113.,
                                123.,123.,123.,123.,123.,123.,123.,
                                133.,133.,133.,133.,133.,133.,133.
}}; // cm

std::vector<int> numberOfTubes{{766, 766, 766, 0, 766, 766, 766,
                                846, 846, 846, 0, 846, 846, 846,
                                927, 927, 927, 0, 927, 927, 927,
                                1007,1007,1007, 0,1007,1007,1007,
                                1087,1087,1087, 0,1087,1087,1087
}};

std::vector<int> layerType{{4, 5, 6, 2, 4, 5, 6,
                            4, 5, 6, 2, 4, 5, 6,
                            4, 5, 6, 2, 4, 5, 6,
                            4, 5, 6, 2, 4, 5, 6,
                            4, 5, 6, 2, 4, 5, 6
}};



double radLenSilicon = 9.370;  // cm
double radLenAluminium = 8.897; // cm

// double thickness{0.00112*radLenSilicon};
// double thickness = 0.02; // cm
//double thickness = 0.01; // cm
double thickness = 0.004; // cm
double rMinStation = 35.7; // cm
double rMaxStation = 135.; // cm

double eps = 1e-10; // cm
double radLenFractionROC = 0.25;
double thicknessROC = radLenFractionROC*radLenAluminium;
double zDrift = 163; // cm
double zROC = zDrift + thicknessROC/2.;
double rMinROC = 27; // cm
double rMaxROC = 141; // cm
 
double radLenFractionFrame = 0.85;
double thicknessFrame = radLenFractionFrame*radLenAluminium;

double zFrame = zDrift + thicknessROC + thicknessFrame/2. + 0.1; // cm
double zFrameCircum1 = zFrame;
double rMinFrameCircum1 = 35; // cm
double rMaxFrameCircum1 = 42; // cm

double zFrameCircum2 = zFrame;
double rMinFrameCircum2 = 120; // cm
double rMaxFrameCircum2 = 141; // cm

double zFrameRadial = zFrame;
double halfYFrameRadial = 4.0;  // cm
double halfXFrameRadial = (sqrt(rMinFrameCircum2*rMinFrameCircum2-halfYFrameRadial*halfYFrameRadial)-rMaxFrameCircum1)/2.-2*eps;
double cXFrameRadial0 = rMaxFrameCircum1 + halfXFrameRadial + eps;
const int nSectors = 12;


double fwdRMin = 0.;
double fwdRMax = rMaxFrameCircum2 + eps;
double fwdHalfZ = positions.back()+10; // cm


#endif