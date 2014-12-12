#!/bin/bash
files=$@

echo "Number of input files: $#"

for f in $files; do
  root=${f%.out}

#-- file exists
  if [ ! -f "$root.wlgrid" ]; then
   echo "--- no such file ---: $root.wlgrid"
   continue
  fi

  mv -i $root.grid $root.grd_grid
  mv -iv $root.wlgrid $root.flx_grid
  mv -i $root.totals $root.tot_energy
  mv -i $root.Lum $root.flx_luminos
  mv -i $root.LumNum $root.flx_lumnum
  mv -i $root.devLum $root.flx_lumdev
  mv -i $root.gamLum $root.flx_gamluminos

  mv -i $root.methodswap $root.grd_methodswap
  mv -i $root.temp $root.grd_temp
  mv -i $root.eraddens $root.grd_eraddens
  mv -i $root.gamdep $root.grd_gamdep
  mv -i $root.capgrey $root.grd_capgrey
  mv -i $root.sig $root.grd_sig

done
