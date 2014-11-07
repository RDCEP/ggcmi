type file;

app (file o) get_inputs () {
   inputs1 stdout = @o;  
}

app detrend(string infile, int latidx, int lonidx, int numlat, int numlon, int nbtlat, int nbtlon, string vrlist, string outdir) {
   detrend infile latidx lonidx numlat numlon nbtlat nbtlon vrlist outdir;
}

type Inputs {
   string infile;
   string vrlist;
   string outdir;
}

int numlat = 360;
int numlon = 720;
int nbtlat = 8;
int nbtlon = 8;

file ff <"finder.out">;
ff = get_inputs();
Inputs io[] = readData(ff);

foreach i in io {
   foreach latidx in [1 : nbtlat] {
      foreach lonidx in [1 : nbtlon] {
         detrend(i.infile, latidx, lonidx, numlat, numlon, nbtlat, nbtlon, i.vrlist, i.outdir);
      }
   }
}
