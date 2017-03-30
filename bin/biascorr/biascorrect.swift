type file;

app (file o) get_inputs (string params) {
   inputs params stdout = @o;
}

app (file o) biascorrect(string inputfile, string reffile, string agglvl, string outdir, string params) {
    biascorrect "-i" inputfile "-r" reffile "-a" agglvl "-o" outdir "-p" params stdout = @o;
}

type Inputs {
    string inputfile;
    string reffile;
    string agglvl;
    string outdir;
}

file ff <"finder.out">;
string params = arg("params");
ff = get_inputs(params);
Inputs irao[] = readData(ff);

foreach i, idx in irao {
    file logfile <single_file_mapper; file = strcat("logs/log_", idx, ".txt")>;
    logfile = biascorrect(i.inputfile, i.reffile, i.agglvl, i.outdir, params);
}
