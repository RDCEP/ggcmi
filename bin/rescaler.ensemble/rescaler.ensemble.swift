type file;

app (file o) get_inputs (string params) {
   inputs params stdout = @o;
}

app (file o) rescaler(string infile, string indir, string mkfile, string agglvl, string outfile, string params) {
   rescaler "-i" infile "-d" indir "-m" mkfile "-a" agglvl "-o" outfile "-p" params stdout = @o;
}

string params = arg("params");

type Inputs {
   string infile;
   string indir;
   string mkfile;
   string agglvl;
   string outfile;
}

file ff <"finder.out">;
ff = get_inputs(params);
Inputs inp[] = readData(ff);

foreach i, idx in inp {
   file logfile <single_file_mapper; file = strcat("logs/log_", idx, ".txt")>;
   logfile = rescaler(i.infile, i.indir, i.mkfile, i.agglvl, i.outfile, params);
}
