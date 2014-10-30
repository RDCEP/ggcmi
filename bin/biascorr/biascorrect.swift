type file;

app (file o) get_inputs () {
   inputs stdout = @o;  
}

app biascorrect(string inputfile, string reffile, string agglvl, string outdir) {
    biascorrect "-i" inputfile "-r" reffile "-a" agglvl "-o" outdir;
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
   file aggfiles[] <filesys_mapper; location = i.indir, pattern = "*">;
   foreach f in aggfiles {
      string fn[] = @strsplit(@f, "__root__");
      biascorrect(fn[1], i.reffile, i.agglvl, i.outdir);
   }
}
