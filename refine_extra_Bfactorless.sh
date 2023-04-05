#!/bin/bash
set -e

#unified script for all

###unique name
ident="1ps_Bfactorless"

###define time points and AF factors to treat. This can produce quite a lot of iterations
declare -a tps=("1ps")

###set dark model to us
dark="dark.pdb"

###set param file to use
param_file="CBD_phenix_no_Bfactor.params"


###wxpert settings
#set dark refine.mtz to use for (calc dmap scrpt)
refine=refine.mtz
echo $refine

declare -a AFs=("25")

#which cif file to use
cif="LBV_modified.cif"
#cif="LBV.cif"


###create all directories in this one
basedir="$PWD"
#basedir=$PWD

###should contain the params file, and the dark files in /dark and extrapolated in extrapol_mtz
sourcedir="$PWD"

###
toggle_real_space_refinement="0" # set to 1 real sapce is not tested yet. The idea is from xplor8 - first reciprocal space to get new phases, compute a 2Fext-Fc map and then real space refine. Good idea, but currently we only refine a part of the molecule in rec space, so not sure if this works great. 

#madan



######### prepare all directories first. This makes sure that all files are present before phenix is run. 

mkdir -p "${basedir}/${ident}/"
cp "$0" "${basedir}/${ident}/" #copy script into output directory
cp $param_file "${basedir}/${ident}/" #also parameter file

for tp in "${tps[@]}" 
do

for AF in "${AFs[@]}"
   
do

  tdir=${basedir}/${ident}/${tp}_$AF  
  echo "in folder" $tdir
  
  #make a folder - check if folder is there, if not make one
  mkdir -p $tdir

  #copy files into folder
  cp $sourcedir/dark/$dark $tdir/dark.pdb
  cp $sourcedir/dark/$cif $tdir
  cp $sourcedir/dark/$refine $tdir/refine.mtz

  cp $sourcedir/extrapolated_mtz/${tp}_extra_${AF}.mtz $tdir
  cp $param_file $tdir/CBD_phenix_no_Bfactor.params
  cp $sourcedir/scripts/mock-dark $tdir

  done
 #make a run file for coot,  one for each timepoint
 coottp=run_coot_${tp}.py
# echo 'set_nomenclature_errors_on_read("auto-correct")' >> $coottp
# echo 'read_cif_dictionary("${tdir}$cif")' >> $coottp 
# echo 'handle_read_draw_molecule("${tdir}dark.pdb")” >> $coottp

done

coot="run_coot_all.py"
#echo 'set_nomenclature_errors_on_read("auto-correct")' >> $coot
#echo 'read_cif_dictionary("${tdir}$cif")' >> $coot
#echo 'handle_read_draw_molecule("${tdir}dark.pdb")” >> $coot

#make directory to collect structures and difference maps
pdir=${basedir}/${ident}/pdbs
mdir=${basedir}/${ident}/maps
ldir=${basedir}/${ident}/logs
mkdir -p $pdir
mkdir -p $mdir
mkdir -p $ldir

#prepare the pcc_all.tzt file
echo "FILE        BproA BproB CproA CproB D207A D207B DA4   DB4   H260A H260B MWA2  PWA3  PWA9  PWB3  PWB9  Rwstar Rfstar Rwfina Rffina " >$pdir/pcc_rec_all.txt
echo "FILE        BproA BproB CproA CproB D207A D207B DA4   DB4   H260A H260B MWA2  PWA3  PWA9  PWB3  PWB9  Rwstar Rfstar Rwfina Rffina " >$pdir/pcc_real_all.txt


#run phenix in target directories. Loop again through the  time points and AFs.
for tp in "${tps[@]}"
do
 for AF in "${AFs[@]}"
 do
  tdir=${basedir}/$ident/${tp}_$AF
  prefix=${ident}_${tp}_${AF}
  cd $tdir
  #phenix.refine extra_$tp_12.mtz dark.pdb CBD_phenix.params
  if [ ! -f "${prefix}_001.pdb" ]; 
   then
   echo "running phenix reciprocal space in $tdir"
   ${seb_phenix}phenix.refine ${tp}_extra_${AF}.mtz monomers.file_name=$cif xray_data.labels="F_${tp},SIGF_${tp}" xray_data.r_free_flags.file_name="${tp}_extra_${AF}.mtz" output.prefix=$prefix CBD_phenix_no_Bfactor.params --overwrite
   else
   echo "not running reciprocal space phenix in $tdir
 ((force run by deleting $tdir${prefix}_001.pdb))"
   fi

  if [ $toggle_real_space_refinement == '1' ];
  then
    pdb_real="${prefix}_001_real_space_refined_000.pdb" 
    if [ ! -f $pdb_real ]; 
     then
      echo "running phenix real-space in $tdir"
      ${seb_phenix}phenix.real_space_refine ${prefix}_001.mtz label='2FOFCWT,PH2FOFCWT' ${prefix}_001.pdb $cif CBD_real_space.params --overwrite   
     else 
     echo "not running real space phenix in $tdir 
 ((force run by deleting $tdir${pdb_real}))"
    fi
  fi

  #make caclulcated difference map
  echo "running calc_dmap_b2.sh in $tdir on ${prefix}_001.pdb \n"
  $sourcedir/scripts/calc_dmap_b2.sh ${prefix}_001.pdb ${tp}_extra_${AF}.mtz ${tp} > $tdir/calc_dmap_b2.log

  echo "running pcc_all.sh in $tdir \n"
  $sourcedir/scripts/pcc_all.sh ${prefix}_001.pdb_wd.map ${sourcedir}/full_dmaps/${tp}_wd.map
  
  rm -rf pcc_out_rec
  mv pcc_out pcc_out_rec  
  mv pcc_summary.txt pcc_rec_summary.txt

  #collect output from  pcc_rec_summary.txt
  printf "%-12s" ${tp}_${AF} >>   $pdir/pcc_rec_all.txt
  grep -o '[0-9].[0-9][0-9][0-9]' pcc_rec_summary.txt | tr '\n' ' ' >> $pdir/pcc_rec_all.txt
  
  #add the Rfree and Rwork from phenix refinement
  tail -3 ${prefix}_001.log | egrep -o '[0-9]{1}.[0-9]{4}'| tr '\n' ' ' >>   $pdir/pcc_rec_all.txt
  echo ' ' >>   $pdir/pcc_rec_all.txt
  ln -f ${prefix}_001.pdb $pdir
  ln -f ${prefix}_001.pdb_wdex.map $mdir
  ln -f ${prefix}_001.log $ldir


  if [ $toggle_real_space_refinement == '1' ];
  then
    echo "running calc_dmap_b2.sh in $tdir on ${pdb_real}"
    $sourcedir/scripts/calc_dmap_b2.sh $pdb_real ${tp}_extra_${AF}.mtz ${tp} > $tdir/calc_dmap_b2_real_space.log
    echo "running pcc_all.sh in $tdir"
    $sourcedir/scripts/pcc_all.sh ${pdb_real}_wd.map ${sourcedir}/full_dmaps/${tp}_wd.map

    rm -rf pcc_out_real
    mv pcc_out pcc_out_real
    mv pcc_summary.txt pcc_real_summary.txt


    #collect output from  pcc_real_summary.txt

    printf "%-12s" ${tp}_${AF} >>   $pdir/pcc_real_all.txt
    grep -o '[0-9].[0-9][0-9][0-9]' pcc_real_summary.txt | tr '\n' ' ' >> $pdir/pcc_real_all.txt
    echo ' ' >>   $pdir/pcc_real_all.txt
    ln -f ${pdb_real} $pdir
    ln -f ${pdb_real}_wdex.map $mdir   
    ln -f ${prefix}_001_real_space_refined_000.log $ldir
  fi

   done
done

echo -en "\007"
echo -en "\007"
echo -en "\007"
echo -en "\007"

#fill the coot command file
#echo "read_map(out_001.pdb)


#fill the coot command file 
