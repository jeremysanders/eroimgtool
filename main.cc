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
#include "event_mode.hh"

int main(int argc, char** argv)
{
  // whether to run in imaging or exposure mode
  std::map<std::string, Pars::runmodetype> modemap{
    {"image", Pars::IMAGE},
    {"expos", Pars::EXPOS},
    {"event", Pars::EVENT}
  };

  // map projection mode names to enum values
  std::map<std::string, Pars::projmodetype> projmodemap{
    {"fov", Pars::AVERAGE_FOV},
    {"fov_sky", Pars::AVERAGE_FOV_SKY},
    {"full", Pars::AVERAGE_FULL},
    {"det", Pars::WHOLE_DET},
    {"radial", Pars::RADIAL},
    {"radial_sym", Pars::RADIAL_SYM},
    {"box", Pars::BOX},
  };

  CLI::App app{"Make eROSITA unvignetted detector exposure maps and images"};
  argv = app.ensure_utf8(argv);

  Pars pars;

  app.add_option("--sources", pars.sources, "List of RA,Dec for sources")
    ->delimiter(',')->expected(1,1000)->required();
  app.add_option("--proj", pars.projmode, "Projection mode")
    ->transform(CLI::CheckedTransformer(projmodemap, CLI::ignore_case))
    ->capture_default_str();
  app.add_option("--proj-args", pars.projargs, "List of arguments for projection")
    ->delimiter(',');
  app.add_option("--tm", pars.tm, "TM number")
    ->check(CLI::Range(1,7))
    ->capture_default_str();
  app.add_option("--pixsize", pars.pixsize, "Pixel size (detector pixels)")
    ->capture_default_str();
  app.add_option("--mask", pars.mask_fn, "Input mask filename")
    ->check(CLI::ExistingFile);
  app.add_option("--mask-pts", pars.maskpts, "Extra masks (list ra,dec,rad_pix)")
    ->delimiter(',');
  app.add_flag("--detmap", pars.detmapmask, "Add CALDB DETMAP mask");
  app.add_flag("--shadowmask", pars.shadowmask, "Add shadow DETMAP mask");
  app.add_option("--gti", pars.gti_fn, "Additional GTI file to merge")
    ->check(CLI::ExistingFile);
  app.add_option("--xw", pars.xw, "X output image size")
    ->capture_default_str();
  app.add_option("--yw", pars.yw, "Y output image size")
    ->capture_default_str();
  app.add_option("--pi-min", pars.pimin, "Minimum PI value (image/event mode)")
    ->capture_default_str();
  app.add_option("--pi-max", pars.pimax, "Maximum PI value (image/event mode)")
    ->capture_default_str();
  app.add_option("--delta-t", pars.deltat, "Time step (s)")
    ->capture_default_str();
  app.add_option("--samples", pars.samples, "Activate sampling for exposure map with number of samples given")
    ->capture_default_str();
  app.add_option("--threads", pars.threads, "Number of threads")
    ->capture_default_str();
  app.add_option("--bitpix", pars.bitpix, "How many bitpix to use for output exposure maps")
    ->capture_default_str();

  app.add_option("mode", pars.mode, "Program mode")
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
      switch(pars.mode)
        {
        case Pars::IMAGE:
          imageMode(pars);
          break;
        case Pars::EXPOS:
          exposMode(pars);
          break;
        case Pars::EVENT:
          eventMode(pars);
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
