type file;

app (file o) get_inputs () {
   inputs stdout = @o;  
}

app rescaler(string infile, string indir, string mkfile, string crmthd, string agglvl, string outfile) {
   rescaler "-i" infile "-d" indir "-m" mkfile "-c" crmthd "-a" agglvl "-o" outfile;
}

type Inputs {
   string infile;
   string indir;
   string mkfile;
   string crmthd;
   string agglvl;
   string outfile;
}

file ff <"finder.out">;
ff = get_inputs();
Inputs inp[] = readData(ff);

foreach i in inp {
   rescaler(i.infile, i.indir, i.mkfile, i.crmthd, i.agglvl, i.outfile);
}
