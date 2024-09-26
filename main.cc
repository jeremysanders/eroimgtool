#include <cmath>
#include <cstdio>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>

#include "CLI11.hpp"

#include "pars.hh"

#include "image_mode.hh"
#include "expos_mode.hh"

int main(int argc, char** argv)
{
  // whether to run in imaging or exposure mode
  enum Mode : int { IMAGE, EXPOS };
  std::map<std::string, Mode> modemap{{"image", IMAGE}, {"expos", EXPOS}};
  Mode mode;

  std::map<std::string, Pars::projmodetype> projmodemap{
    {"fov", Pars::AVERAGE_FOV},
    {"fov_sky", Pars::AVERAGE_FOV_SKY},
    {"det", Pars::WHOLE_DET}};

  CLI::App app{"Make eROSITA unvignetted detector exposure maps"};
  argv = app.ensure_utf8(argv);

  Pars pars;
  app.add_option("--proj", pars.projmode, "Projection mode")
    ->transform(CLI::CheckedTransformer(projmodemap, CLI::ignore_case))
    ->capture_default_str();
  app.add_option("--proj-args", pars.projargs, "List of arguments for projection");
  app.add_option("--tm", pars.tm, "TM number")
    ->check(CLI::Range(1,7))
    ->capture_default_str();
  app.add_option("--ra", pars.src_ra, "RA of source (deg)")
    ->required();
  app.add_option("--dec", pars.src_dec, "Dec of source (deg)")
    ->required();
  app.add_option("--pixsize", pars.pixsize, "Pixel size (detector pixels)")
    ->capture_default_str();
  app.add_option("--mask", pars.mask_fn, "Input mask filename")
    ->check(CLI::ExistingFile);
  app.add_option("--xw", pars.xw, "X output image size")
    ->capture_default_str();
  app.add_option("--yw", pars.yw, "Y output image size")
    ->capture_default_str();
  app.add_option("--pi-min", pars.pimin, "Minimum PI value (image mode)")
    ->capture_default_str();
  app.add_option("--pi-max", pars.pimax, "Maximum PI value (image mode)")
    ->capture_default_str();
  app.add_option("--delta-t", pars.deltat, "Time step (s)")
    ->capture_default_str();
  app.add_option("--threads", pars.threads, "Number of threads")
    ->capture_default_str();

  app.add_option("mode", mode, "Program mode")
    ->required()
    ->transform(CLI::CheckedTransformer(modemap, CLI::ignore_case));
  app.add_option("event", pars.evt_fn, "Event filename")
    ->required()
    ->check(CLI::ExistingFile);
  app.add_option("image", pars.out_fn, "Output image filename")
    ->required();

  CLI11_PARSE(app, argc, argv);

  // pars.tm = 3;
  // pars.evt_fn = "em01_258138_020_ML00001_004_c946/evt.fits.gz";
  // pars.mask_fn = "em01_258138_020_ML00001_004_c946/030_mask_final.fits.gz";
  // pars.out_fn = "test.fits";
  // pars.src_ra = 57.3466206;
  // pars.src_dec = -11.9909090;
  // pars.pixsize = 1/8.f;
  // pars.threads = 4;
  // pars.src_ra = 255.7054905; //57.3466206;
  // pars.src_dec = -48.7900230; //-11.9909090;

  // pars.xw = 1024;
  // pars.yw = 1024;

  // pars.xw *= 4;
  // pars.yw *= 4;
  // pars.pixsize /= 4;

  //pars.deltat = 0.1;
  //pars.projmode = Pars::WHOLE_DET;
  //pars.projmode = Pars::AVERAGE_FOV_SKY;

  imageMode(pars);
  //exposMode(pars);

  return 0;
}
