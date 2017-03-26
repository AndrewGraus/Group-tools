This folder contains programs to modify the gizmo hdf5 outputs so you can actually run halo finders on them
You need to do this no matter which halo finder you run because

1) For AHF it currently only runs on gadget binaries, so you need to convert the hdf5 files back to binaries

2) For rockstar you need to add the hi res particle mass to the header or rockstar will think the particles are
  massless
  
  NOTE: AHF seems to be giving bizzare results, or just not running at all, at times, and then at other times runs
  just fine.  We have not tracked down why this is so use AHF at your own risk, and DEFINATELY verify anything with
  rockstar.  I'm going to be honest, you should run both always reguardless
