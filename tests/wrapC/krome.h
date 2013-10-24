#ifndef KROME_H_
#define KROME_H_

extern int krome_idx_H = 0; //H
extern int krome_idx_E = 1; //E
extern int krome_idx_Hj = 2; //H+
extern int krome_idx_HE = 3; //HE
extern int krome_idx_HEj = 4; //HE+
extern int krome_idx_HEjj = 5; //HE++
extern int krome_idx_Hk = 6; //H-
extern int krome_idx_H2 = 7; //H2
extern int krome_idx_H2j = 8; //H2+
extern int krome_idx_HEHj = 9; //HEH+
extern int krome_idx_H3j = 10; //H3+
extern int krome_idx_CR = 11; //CR
extern int krome_idx_g = 12; //g
extern int krome_idx_Tgas = 13; //Tgas
extern int krome_idx_dummy = 14; //dummy

extern char* krome_names[] = {
  "H",
  "E",
  "H+",
  "HE",
  "HE+",
  "HE++",
  "H-",
  "H2",
  "H2+",
  "HEH+",
  "H3+",
  "CR",
  "g",
  "Tgas",
  "dummy"
};
extern int krome_nmols = 11; //number of species
extern int krome_nrea = 36; //number of reactions

#endif
