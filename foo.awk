#!/usr/local/bin/awk -f

{
    ## FNR is the line number of the current file, NR is the number of 
    ## lines that have been processed. If you only give one file to
    ## awk, FNR will always equal NR. If you give more than one file,
    ## FNR will go back to 1 when the next file is reached but NR
    ## will continue incrementing. Therefore, NR == FNR only while
    ## the first file is being processed.
    if(NR == FNR){
      ## If this is the first file, save the values of $1
      ## in the array n.
      n[$1] = 0
    }
    ## If we have moved on to the 2nd file
    else{
      ## If the 3rd field of the second file exists in
      ## the first file.
      if($3 in n){
        ## Add the value of the 5th field to the corresponding value
        ## of the n array.
        n[$3]+=$5
      }
    }
}
## The END{} block is executed after all files have been processed.
## This is useful since you may have more than one line whose 3rd
## field was specified in the first file so you don't want to print
## as you process the files.
END{
    ## For each element in the n array
    for (i in n){
    ## print the element itself and then its value
    print i,":",n[i];
    }
}