input=$1
output=$2

if [ -f $output ]; then \rm $output; fi

if [ -r $input ];
then echo $input >> $output
fi
