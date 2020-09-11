#include <ctime>
#include <string>

#include "imufile.h"
#include "ins.h"
#include "kf.h"

using nstd::navigation::ImuErr;
using nstd::navigation::ImuFile;
using nstd::navigation::Ins;
using nstd::navigation::Kf;
using nstd::navigation::Pva;

int main(int argc, char const* argv[]) {
  std::string inpath(
      "/home/long/Documents/intnavdata/stim300/STIM30020161026_03.bin");
  std::string gnss_input_path(
      "/home/long/Documents/intnavdata/stim300/trueob.txt");
  // "/home/long/Documents/intnavdata/stim300/GBRTK_edit.txt");
  constexpr double pos[3]{30.4052867239, 114.2679271285, 24};
  constexpr double fs = 125;
  constexpr double gyrobias[3]{0.5, 0.5, 0.5};
  constexpr double accebias[3]{0.05 * GRAVITY * 0.001, 0.05 * GRAVITY * 0.001,
                               0.05 * GRAVITY * 0.001};
  constexpr double ARW[3]{0.15 * DEG_TO_RAD / 60, 0.15 * DEG_TO_RAD / 60,
                          0.15 * DEG_TO_RAD / 60};
  constexpr double VRW[3]{0.06 / 60, 0.06 / 60, 0.06 / 60};
  double obvar[6]{4 / RE, 4 / RE, 4, 0.1, 0.1, 0.1};
  double att[3]{0};
  double pos_cp_rad[3]{0};
  cblas_dcopy(3, pos, 1, pos_cp_rad, 1);
  cblas_dscal(3, DEG_TO_RAD, pos_cp_rad, 1);
  ImuFile data(inpath);
  nstd::navigation::Alignment(data, pos, &att);
  double vel[3]{0, 0, 0};
  ImuErr err(gyrobias, accebias, ARW, VRW);
  Pva inipva(pos, vel, att);
  Ins inertial(inipva, fs, err);
  Kf kalman(err, obvar);
  nstd::navigation::Gnss sat;
  sat.ReadFile(gnss_input_path);
  // kalman.SetGnssBreakTime(std::vector<double>({287200, 287260}));
  kalman.Update(inertial, sat, data);
  std::cout << "saving file" << std::endl;
  data.Save("resource/test.bin");
  nstd::matrix::SaveMat("resource/state.txt", kalman.x_total_.data(), 10,
                        kalman.x_total_.size() / 10);
  return 0;
}