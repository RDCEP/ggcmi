type file;

app (file o) get_inputs () {
   inputs2 stdout = @o;  
}

app combine(string indir, int latidx, string infile) {
   combine indir latidx infile;
}

type Inputs {
   string indir;
   string infile;
}

int nbtlat = 8;

file ff <"finder.out">;
ff = get_inputs();
Inputs ii[] = readData(ff);

foreach i in ii {
   foreach latidx in [1 : nbtlat] {
      combine(i.indir, latidx, i.infile);
   }
}
