#pragma once

#include <string>

#include "global.h"
#include "hdf5.h"

class Reader
{
private:
  void read_face_ctr(hid_t fid);
  void read_face_normal(hid_t fid);
  void read_face_area(hid_t fid);
  void read_data(hid_t fid);
  void read_input();

  string input_fnameS;
  string fwh_surf_fname;

public:
  Reader(char *input_fnameC);
  void read();

  ~Reader();
};
