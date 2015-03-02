#!/bin/sh

if ! [ -z $1 ]
then
  echo "#!/bin/sh" >> .submission
  echo "#" >> .submission
  echo "#PBS -N arj32"  >> .submission
  if [ -z $2 ]
  then
    echo "#PBS -l nodes=1" >> .submission
  else
    echo "#PBS -l nodes=compute-0-"$2"" >> .submission
  fi
  echo  >> .submission
  echo  >> .submission
  echo "module add R/3.0.1" >> .submission
  echo "cd /share/scratch/arj32/tcga-analysis" >> .submission
  echo "R CMD BATCH" $1 >> .submission

  qsub .submission
  rm .submission
fi
