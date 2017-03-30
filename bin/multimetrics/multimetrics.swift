type file;

app (file o) get_inputs (string params) {
   inputs params stdout = @o;
}

app (file o) multimetrics(string inputfile, string reffile, string agglvl, string outdir, string params) {
    multimetrics inputfile reffile agglvl outdir params stdout = @o;
}

type Inputs {
    string inputfile;
    string reffile;
    string agglvl;
    string outdir;
}

string params = arg("params");
file ff <"finder.out">;
ff = get_inputs(params);
Inputs irao[] = readData(ff);

foreach i, idx in irao {
   file logfile <single_file_mapper; file = strcat("logs/log_", idx, ".txt")>;
   logfile = multimetrics(i.inputfile, i.reffile, i.agglvl, i.outdir, params);
}
