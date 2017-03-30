type file;

app (file o) get_inputs (string params) {
   inputs params stdout = @o;
}

app modelensemble(string indir, string metricsdir, string agglvl, string weather, string crop, string metric, string outdir) {
    modelensemble "-d" indir "-m" metricsdir "-a" agglvl "-w" weather "-c" crop "--metric" metric "-o" outdir;
}

type Inputs {
    string indir;
    string metricsdir;
    string agglvl;
    string outdir;
}

string climates[] = strsplit(arg("w"), ",");
string crops[]    = strsplit(arg("c"), ",");
string metrics[]  = strsplit(arg("m"), ",");
string params = arg("params");
file ff <"finder.out">;
ff = get_inputs(params);
Inputs inp[] = readData(ff);

foreach i in inp {
   foreach w in climates {
      foreach c in crops {
         foreach m in metrics {
            modelensemble(i.indir, i.metricsdir, i.agglvl, w, c, m, i.outdir);
         }
      }
   }
}
