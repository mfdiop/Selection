#!/bin/sh

#SBATCH -A mdiop
#SBATCH -J filter
#SBATCH -p bigmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1   
#SBATCH --mem=50000
#SBATCH --array=1-10
#SBATCH --error="filter.error"
#SBATCH --output="filter-%A_%a.out"    
#SBATCH --mail-type=BEGIN,FAIL,END          
#SBATCH --mail-user=mdiop@mrc.gm

typeset -F SECONDS
module load R

#Treat expanding unset parameters as error 
set -u
#Exit with status 1 if a command throws an error
set -e

if [ $# -le 0 ]
then
  echo ""
  echo "Filter missing data from vcf file"
  echo "Dependencies: vcf files / working directory"
  echo "Produces: filtered.vcf.gz"
  echo ""
  echo "usage: PATH/scripts/filter.sh -p path -m maf"
  echo ""
  echo "-p, --path:  Full path to folder where VCF files are stored." 
  echo "-m, --maf:   Minor Allele Frequency."
  echo ""
  exit 1;
fi

PATH_TO_SCRIPT=$0
PATH_TO_SCRIPT=$(echo ${PATH_TO_SCRIPT} | awk -F\scripts/filter.sh '{print $1}')

##############################################################

########### Read arguments from command line #################

POSITIONAL=()
while [[ $# -gt 0 ]]
do
  key="$1"

  case $key in
    -p|--path)
      path="$2"
      shift # past argument
      shift # past value
      ;;
	-m|--maf)
      maf="$2"
      shift # past argument
      shift # past value
      ;;
    *)    # unknown option
      POSITIONAL+=("$1") # save it in an array for later
      shift # past argument
      ;;
  esac
done


echo "********************************"
echo "Parameters passed to script:"
echo "Full Path           = $path"
echo "Minor Allele Freq   = $maf"
echo "********************************"

#### check if data exists
check_file_existence (){
  if [ ! -f $1 ]; then
    echo "0"
  else
    echo "1"
  fi
}

#####################################
vcf_files=`ls $path/*.vcf.gz`
vcf=`echo $vcf_files | cut -f $SLURM_ARRAY_TASK_ID -d" "`

Rscript filter.R $path $vcf $maf

sleep 1
echo -e "\n\nWall time is "$SECONDS" seconds"













