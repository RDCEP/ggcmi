type file;

app (file o) get_inputs () {
   inputs stdout = @o;  
}

app multiensemble(string inputfile, string reffile, string agglvl, string outdir) {
    multiensemble inputfile reffile agglvl outdir;
}

type Inputs {
    string indir;
    string reffile;
    string agglvl;
    string outdir;
}

file ff <"finder.out">;
ff = get_inputs();
Inputs iro[] = readData(ff);

foreach i in iro {
   file esfiles[] <filesys_mapper; location = i.indir, pattern = "*">;
   foreach f in esfiles {
      string fn[] = strsplit(@f, "file://localhost/");
      multiensemble(fn[1], i.reffile, i.agglvl, i.outdir);
   }
}
