type file;

app (file o) get_inputs () {
   inputs stdout = @o;  
}

app (file o) downscaler(string irfile, string rffile, string agfile, string reffile, string mkfile, string wtfile, string agglvl, string vname, string outfile) {
   downscaler "-i" irfile "-r" rffile "-a" agfile "-f" reffile "-m" mkfile "-w" wtfile "--agglvl" agglvl "-v" vname "-o" outfile stdout = @o;
}

type Inputs {
   string irfile;
   string rffile;
   string agfile;
   string wtfile;
   string vname;
   string outfile;
}

string mkfile  = arg("mkfile");
string reffile = arg("reffile");
string agglvl  = arg("agglvl");

file ff <"finder.out">;
ff = get_inputs();
Inputs inp[] = readData(ff);

foreach i, idx in inp {
   file logfile <single_file_mapper; file = strcat("logs/log_", idx, ".txt")>;
   logfile = downscaler(i.irfile, i.rffile, i.agfile, reffile, mkfile, i.wtfile, agglvl, i.vname, i.outfile);
}
