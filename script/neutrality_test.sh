#!/bin/sh

#SBATCH -A mdiop
#SBATCH -J run
#SBATCH -p bigmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1   
#SBATCH --mem=5000
#SBATCH --error="run.error"
#SBATCH --output="run.out"    
#SBATCH --mail-type=BEGIN,FAIL,END          
#SBATCH --mail-user=mdiop@mrc.gm

typeset -F SECONDS
module load R
module load bcftools
module load tabix

#Treat expanding unset parameters as error 
set -u
#Exit with status 1 if a command throws an error
set -e

if [ $# -le 0 ]
then
  echo ""
  echo "Compute neutrality test from the filtered vcf file"
  echo "Dependencies: filtered vcf files / gff file/ Output directory"
  echo "Produces: Tajima.txt"
  echo ""
  echo "usage: PATH/scripts/neutrality_test.sh -p path -g gff -o output"
  echo ""
  echo "-p, --path:    Full path to folder where filtered VCF files are stored."
  echo "-g, --gff:     Full path to Plasmodium falciparum gff file." 
  echo "-g, --output:  Full path to output folder." 
  echo ""
  exit 1;
fi

PATH_TO_SCRIPT=$0
PATH_TO_SCRIPT=$(echo ${PATH_TO_SCRIPT} | awk -F\scripts/neutrality_test.sh '{print $1}')

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
	-g|--gff)
      gff="$2"
      shift # past argument
      shift # past value
      ;;
	-o|--output)
      output="$2"
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
echo "Directory           = $path"
echo "GFF file            = $gff"
echo "Output              = $output"
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

# Load VCF files
vcf_files=`ls $path/*.maf4.vcf`

# Get single VCF for parallel computing
vcf=`echo $vcf_files | cut -f $SLURM_ARRAY_TASK_ID -d" "`

# Run R script
Rscript neutrality_test.R $vcf $gff $output

sleep 1
echo -e "\n\nWall time is "$SECONDS" seconds"














