type file;

app (file o) get_inputs () {
   inputs stdout = @o;  
}

app (file o) biascorrect(string inputfile, string reffile, string agglvl, string outdir) {
    biascorrect "-i" inputfile "-r" reffile "-a" agglvl "-o" outdir stdout = @o;
}

type Inputs {
    string inputfile;
    string reffile;
    string agglvl;
    string outdir;
}

file ff <"finder.out">;
ff = get_inputs();
Inputs irao[] = readData(ff);

foreach i, idx in irao {
    file logfile <single_file_mapper; file = strcat("logs/log_", idx, ".txt")>;
    logfile = biascorrect(i.inputfile, i.reffile, i.agglvl, i.outdir);
}
