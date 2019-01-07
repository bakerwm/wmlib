# calculate the base content for fastx file

{
    if(NR%4==2) {
        lenmax=length($1)
        if (max == 0) { max=lenmax }
        for (i=1; i<=max; i++) {
            q=substr($1,i,1)
            s[i][q]++
            #print c
        }
    }
} 
END {
    print "position","A","C","G","T","N"
    for (i in s) {
        if (isarray(s[i])) {
            if (s[i]["A"] == "") {s[i]["A"]=0}
            if (s[i]["C"] == "") {s[i]["C"]=0}
            if (s[i]["G"] == "") {s[i]["G"]=0}
            if (s[i]["T"] == "") {s[i]["T"]=0}
            if (s[i]["N"] == "") {s[i]["N"]=0}
            print i,s[i]["A"],s[i]["C"],s[i]["G"],s[i]["T"],s[i]["N"]
        }
    }
}