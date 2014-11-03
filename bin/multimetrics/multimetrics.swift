type file;

app (file o) get_inputs () {
   inputs stdout = @o;  
}

app multimetrics(string inputfile, string reffile, string agglvl, string outdir) {
    multimetrics inputfile reffile agglvl outdir;
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
   file bcfiles[] <filesys_mapper; location = i.indir, pattern = "*">;
   foreach f in bcfiles {
      string fn[] = @strsplit(@f, "__root__");
      multimetrics(fn[1], i.reffile, i.agglvl, i.outdir);
   }
}
