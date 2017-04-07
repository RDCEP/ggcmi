type file;

app (file o) get_inputs (string params) {
   inputs params stdout = @o;
}

app (file o) rescaler(string params, string irfile, string rffile, string bcfile, string mkfile, string wtfile, string agglvl, string vname, string outfile) {
   rescaler "-p" params "-i" irfile "-r" rffile "-b" bcfile "-m" mkfile "-w" wtfile "-a" agglvl "-v" vname "-o" outfile stdout = @o;
}

type Inputs {
   string irfile;
   string rffile;
   string bcfile;
   string wtfile;
   string vname;
   string outfile;
}

string mkfile = arg("mkfile");
string agglvl = arg("agglvl");
string params = arg("params");

file ff <"finder.out">;
ff = get_inputs(params);
Inputs inp[] = readData(ff);

foreach i, idx in inp {
   file logfile <single_file_mapper; file = strcat("logs/log_", idx, ".txt")>;
   logfile = rescaler(params, i.irfile, i.rffile, i.bcfile, mkfile, i.wtfile, agglvl, i.vname, i.outfile);
}
