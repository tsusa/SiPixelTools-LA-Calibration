#ifndef LAStructures_h
#define LAStructures_h

static const int maxpix = 1000;
struct Pixinfo {
  int npix;
  float row[maxpix];
  float col[maxpix];
  float adc[maxpix];
  float x[maxpix];
  float y[maxpix];
};

struct ColInfo {
  Int_t ncol;
  Int_t dcol[maxpix];
  Float_t adc[maxpix];
  Float_t depth[maxpix];
};

struct Hit {
  float x;
  float y;
  double alpha;
  double beta;
  double gamma;
};

struct Clust {
  float x;
  float y;
  float charge;
  int size_x;
  int size_y;
  int maxPixelCol;
  int maxPixelRow;
  int minPixelCol;
  int minPixelRow;
};
struct Rechit {
  float x;
  float y;
};
struct LAResultsGA {
  float LA_lin;
  float LA_linErr;
  float LA_width;
  float LAwidthErr;
  float LA_width_zeroSub;
  float LA_widthErr_zeroSub;
  float LA_der;
  float LA_derErr; 
  char  name[256];
};

#endif
