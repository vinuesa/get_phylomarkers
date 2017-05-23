#!/usr/bin/env bash

# AUTHOR: Pablo Vinuesa; http://www.ccg.unam.mx/~vinuesa/
# AIM: run commands in parallel to process multiple files using pexec

progname=$(basename $0)
VERSION='v0.2_22May17' # removed the echo $file ... after command for better portability! $command; echo $file ...
            # v0.1 Sept 23rd, 2013


function check_dependencies()
{
    for programname in pexec
    do

       #if which $programname >/dev/null; then <== avoid which
       # see: http://stackoverflow.com/questions/592620/check-if-a-program-exists-from-a-bash-script

       bin=$(type -P $programname)
       if [ -z $bin ]; then
          echo
          echo "# ERROR: $programname not in place!"
          echo "# ... you will need to install \"$programname\" first or include it in \$PATH"
          echo "# ... exiting"
          exit 1
       fi
    done

}



function print_help()
{
    cat << HELP
  
    $progname v $VERSION usage: requires 2 or max 3 args
    $progname <file extension name> <'command'> [no_of_cores]

    example1: $progname faa   'muscle < \$file > \${file%.faa}_musAln.FAA'
    example2: $progname faa   'mafft-linsi \$file > \${file%.faa}_mlinsi.FAA'
    example3: $progname FAA   'FastTree -slow -slownni -spr 25 -sprlength 12 -log \${file%.FAA}_FT.log < \$file > \${file%.FAA}_FT.ph'
    example4: $progname FASTA 'FastTree -nt -gtr -gamma -slow -slownni -spr 25 -sprlength 12 -log \${file%.FASTA}_FT.log < \$file > \${file%.FASTA}_FT.ph'

HELP

    exit 2
}



if [ $# -lt 2 -o $# -gt 3 ]
then
    print_help
fi


check_dependencies

ext=$1
command=$2
no_of_cores=$3


total_files=$(ls *$ext | wc -l)

if [ -z $no_of_cores ]
then
    pexec -r *.$ext -e file -c -o - -- "$command"
else
    pexec -n $no_of_cores -r *.$ext -e file -c -o - -- "$command"
fi

echo 
echo "done processing $total_files files ..."
