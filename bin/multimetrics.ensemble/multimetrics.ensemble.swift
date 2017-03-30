type file;

app (file o) get_inputs (string params) {
   inputs params stdout = @o;
}

app multiensemble(string inputfile, string reffile, string agglvl, string outdir, string params) {
    multiensemble inputfile reffile agglvl outdir params;
}

type Inputs {
    string indir;
    string reffile;
    string agglvl;
    string outdir;
}

string params = arg("params");
file ff <"finder.out">;
ff = get_inputs(params);
Inputs iro[] = readData(ff);

foreach i in iro {
   file esfiles[] <filesys_mapper; location = i.indir, pattern = "*">;
   foreach f in esfiles {
      string fn[] = strsplit(@f, "file://localhost/");
      multiensemble(fn[1], i.reffile, i.agglvl, i.outdir, params);
   }
}
