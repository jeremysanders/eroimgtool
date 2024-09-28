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

  // map projection mode names to enum values
  std::map<std::string, Pars::projmodetype> projmodemap{
    {"fov", Pars::AVERAGE_FOV},
    {"fov_sky", Pars::AVERAGE_FOV_SKY},
    {"det", Pars::WHOLE_DET},
    {"radial", Pars::RADIAL},
    {"radial_sym", Pars::RADIAL_SYM},
  };

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
  app.add_option("--mask-src-rad", pars.masksrcrad, "Source mask radius (pix)")
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

  try
    {
      switch(mode)
        {
        case IMAGE:
          imageMode(pars);
          break;
        case EXPOS:
          exposMode(pars);
          break;
        }
    }
  catch(std::runtime_error& e)
    {
      std::fprintf(stderr, "\nError: %s\n", e.what());
      return 1;
    }

  return 0;
}
